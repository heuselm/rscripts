## Process MaxQuant results
devtools::install_github("bartongroup/proteusLabelFree")
devtools::install_github("bartongroup/Proteus",
build_opts= c("--no-resave-data", "--no-manual"), build_vignettes=FALSE)

# load packages
library(data.table)
library(ggplot2)
library(proteus)
library(shiny)
library(plotly)
library(htmlwidgets)
library(ggrepel)

## Define inputs
# evidence --> quantdata
evidenceFile = "combined/txt/evidence.txt"
ev = fread(evidenceFile)

# metadata
# example
fread(system.file("extdata", "metadata.txt", package="proteusLabelFree"))
# write template
fwrite(data.table("experiment" = unique(ev$Experiment),
                  "measure" = rep("Intensity", length(unique(ev$Experiment))),
                  "sample" = unique(ev$Experiment),
                  "condition" = rep("FILL-IN", length(unique(ev$Experiment))),
                  "replicate" = rep("FILL-IN", length(unique(ev$Experiment)))),
       file = "experimentalDesignMetadataFile_template.txt", sep = "\t")

metadataFile = "experimentalDesignMetadataFile.txt"
meta <- read.delim(metadataFile, header=TRUE, sep="\t")

# Simplify Protein group identifiers
simplification_table = fread("../210622_IEX-hlp-test_fragpipe_analysis/simplification_table.csv")
simplifyCloneIDs = function(cloneIDs_to_simplify, simplification_table){
  cloneIDs_simplified = copy(cloneIDs_to_simplify)
  # subset simplification table to those truly present in cloneIDs_to_simplify
  simplification_table[, present:=length(grep(gsub("\\|","\\\\|", cloneID_old), cloneIDs_to_simplify)), cloneID_old]
  # print(paste("Replacement operations:",nrow(simplification_table[present>0])))
  for (i in 1:nrow(simplification_table[present>0])){
    cloneIDs_simplified = gsub(gsub("\\|","\\\\|", simplification_table[present>0]$cloneID_old[i]),
                               simplification_table[present>0]$cloneID_new[i], cloneIDs_simplified)
  }
  return(data.table(cloneIDs_to_simplify, cloneIDs_simplified))
}
protein_groups_detected = unique(ev$Proteins)

pgsimp = simplifyCloneIDs(protein_groups_detected, simplification_table)

ev = merge(ev, pgsimp, by.x = "Proteins", by.y = "cloneIDs_to_simplify", all.y = F)
ev$Proteins = NULL
setnames(ev, "cloneIDs_simplified", "Proteins")
ev[, `Leading proteins`:=Proteins]
ev[, `Leading razor protein`:=Proteins]

write.table(ev, file = "evidence_simplified.txt", sep = "\t", row.names = F, quote = F)

# Plot IDs
nruns = length(unique(ev$`Raw file`))

ev = merge(ev, meta, by.x = "Experiment", by.y = "experiment")

ids.prot = ev[, length(unique(Proteins)), .(`Raw file`,sample,condition,replicate)]
ids.pep = ev[, length(unique(`Modified sequence`)), .(`Raw file`,sample,condition,replicate)]

ids = rbind(ids.prot[, level:="Protein.groups"],
            ids.pep[, level:="Modified.peptides"])

ggplot(ids, aes(x = condition, y = V1, fill = condition)) + geom_bar(stat = "identity") + ggtitle("Identifications MaxQuant") +
  theme(axis.text.x = element_text(angle = 90)) +
  geom_text(aes(label=V1), vjust = -0.2)+
  ylab("N") +
  facet_wrap(~level, scales = "free_y")
ggsave("01_IDs_manual.pdf", width = (nruns/2)+8)


# Proteus differential analysis
################################
evi = readEvidenceFile(file = "evidence_simplified.txt")

pepdat <- makePeptideTable(evi, meta)
protdat <- makeProteinTable(pepdat)

# Plot ID numbers
plotCount(pepdat)
ggsave("02_pepCount_proteus.pdf")

plotCount(protdat)
ggsave("02_protCount_proteus_min2pep.pdf")

# Plot some example peptides/proteins
pdf("06_protein_plots.pdf", width = nruns/3)
for (i in 1:length(unique(ev$Proteins))){
  pi = plotProtPeptides(pepdat = pepdat, protein = sample(pepdat$proteins, 1))
  plot(pi)
}
dev.off()

# Plot heatmap of protein intensities as summarized by Proteus
library(pheatmap)
protdat$tab[is.na(protdat$tab)] = 0
protdat$tab[is.infinite(protdat$tab)] = 0

matrix = protdat$tab[rowSums(protdat$tab)>0,]

pheatmap(matrix, scale = "row",
         cluster_cols = F,
         main = "Heatmap all runs rowscaled, MaxQuant",
         color = viridis::inferno(30))
dev.copy(pdf, file = "Heatmap_proteus.pdf", width = 18, height = 7)
dev.off()

# write out protein intensity long format
# Write out long table
ev[, Intensity.prot:=sum(Intensity), .(Proteins, Experiment)]
ev$replicate = NULL
ev[, replicate:=gsub("iex-hlp-test-", "", unlist(strsplit(Experiment, split = "_"))[1]), Experiment]
ev[, sample:=unlist(strsplit(Experiment, split = "_"))[2], Experiment]

unique(ev$condition)
unique(ev$sample)
unique(ev$replicate)

data_long = unique(ev[, .(Proteins, replicate, sample, Intensity.prot)])
names(data_long) = c("protein_id", "replicate", "sample", "intensity")
fwrite(data_long, file = "prot_long_intsum.csv")

data_long[, analysis := "MaxQuant"]
data_long_fragpipe = fread("../210622_IEX-hlp-test_fragpipe_analysis/prot_long_intsum.csv")

# Combine datasets
data_long_fragpipe[, analysis:= "FragPipe"]
data_long = rbind(data_long, data_long_fragpipe)
data_long[, sample:=as.numeric(gsub("S","",sample)),sample]

ggplot(data_long[, sum(intensity), .(sample, replicate,analysis)]) +
  geom_bar(aes(x = sample, y = V1, fill = paste(analysis)),
  stat = "identity", position = position_dodge2(preserve = "single", width = 0.5)) +
  facet_wrap(~replicate) + ylab("Total Intensity")

# Draw all lines across both replicates
data_long[, commercial_ab:=FALSE]
data_long[grep("_H$|_L$", protein_id), commercial_ab:=TRUE]

data_long[, sample:=factor(sample, levels = seq(1,48))]

p = ggplot(data_long) +
  geom_line(aes(x = sample, y = intensity, color = protein_id,
                color = protein_id, group = paste(protein_id,replicate,analysis))) +
  geom_point(aes(x = sample, y = intensity, color = protein_id,
                 color = protein_id)) +
  facet_wrap(~analysis+replicate, scales = "free_y") +
  theme(legend.position = "none",
                                          axis.text.x = element_text(angle = 90)) +
  ggtitle("IEX hlp test - All protein profiles")
p
ggplotly(p)

p1 = ggplot(data_long[commercial_ab == TRUE]) +
  geom_line(aes(x = sample, y = intensity, color = protein_id,
                color = protein_id, group = paste(protein_id,replicate,analysis))) +
  geom_point(aes(x = sample, y = intensity, color = protein_id,
                 color = protein_id)) +
  facet_wrap(~analysis+replicate, scales = "free_y") +
  theme(legend.position = "bottom", axis.text.x = element_text(angle = 90)) +
  ggtitle("IEX hlp test - Commercial AB protein profiles")
p1
ggplotly(p1)

p2 = ggplot(data_long[commercial_ab == TRUE]) +
  geom_line(aes(x = sample, y = intensity, color = protein_id,
                color = protein_id, group = paste(protein_id,replicate,analysis))) +
  geom_point(aes(x = sample, y = intensity, color = protein_id,
                 color = protein_id)) +
  facet_wrap(~replicate+analysis, scales = "free_y") +
  theme(legend.position = "bottom", axis.text.x = element_text(angle = 90)) +
  ggtitle("IEX hlp test - Commercial AB protein profiles")
p2
ggplotly(p2)

mat = dcast(data_long, protein_id~analysis+replicate+sample, value.var = "intensity")

mat.m = as.matrix(mat[, 2:ncol(mat)])
row.names(mat.m) = mat$protein_id

mat.m[is.na(mat.m)] = 0
mat.m.log10 = log10(mat.m)
mat.m.log10[is.infinite(mat.m.log10)] = 0

pheatmap(mat.m.log10, show_rownames = F,
        # scale = "row",
         cluster_cols = F,
         gaps_col = c(32,47, 82),
         main = "Heatmap all runs, log10, FragPipe + MaxQuant",
         color = viridis::inferno(30))
dev.copy(pdf, file = "Heatmap_global_log10.pdf", width = 24, height = 14)
dev.off()


