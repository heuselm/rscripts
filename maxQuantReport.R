# maxQuantReport
# mhe 22.06.2021

# Parameters
# -------------------------------------------------------
evidenceFile <- choose.files(caption = "Choose evidence.txt file from MaxQuant output combined/txt/ folder")
quantile_normalization = TRUE

# --------------------------------------------------------

# Dependencies
# devtools::install_github("bartongroup/proteusLabelFree")
# devtools::install_github("bartongroup/Proteus",
# build_opts= c("--no-resave-data", "--no-manual"), build_vignettes=FALSE)

# load packages
library(data.table)
library(ggplot2)
library(proteus)
library(pheatmap)
library(rmarkdown)
library(psych)

## Processing
# evidence.txt
message(evidenceFile)
message(dirname(evidenceFile))
message(basename(evidenceFile))
Sys.sleep(2)

setwd(dirname(evidenceFile))
ev = fread(evidenceFile)

# metadata
# write template file for user
fwrite(data.table("experiment" = unique(ev$Experiment),
                  "measure" = rep("Intensity", length(unique(ev$Experiment))),
                  "sample" = unique(ev$Experiment),
                  "condition" = rep("FILL-IN", length(unique(ev$Experiment))),
                  "replicate" = rep("FILL-IN", length(unique(ev$Experiment)))),
       file = "experimentalDesignMetadataFile_template.csv")

message("experimentalDesignMetadataFile_template.csv written to txt folder - fill in experimental metadata and save as experimentalDesignMetadataFile.csv")

# read filled template back from user
metadataFile = choose.files(caption = "Choose the completed experimentalDesignMetadataFile.tsv file from MaxQuant output combined/txt/ folder")
meta <- fread(metadataFile, header=TRUE)
meta

# Plot IDs
############
nruns = length(unique(ev$`Raw file`))
ev = merge(ev, meta, by.x = "Experiment", by.y = "experiment")

ids.prot = ev[, length(unique(Proteins)), .(`Raw file`,sample,condition,replicate)]
ids.pep = ev[, length(unique(`Modified sequence`)), .(`Raw file`,sample,condition,replicate)]

ids = rbind(ids.prot[, level:="Protein.groups"],
            ids.pep[, level:="Modified.peptides"])

ggplot(ids, aes(x = paste(condition,"R", replicate), y = V1, fill = condition)) + geom_bar(stat = "identity") + ggtitle("Identifications MaxQuant") +
  theme(axis.text.x = element_text(angle = 90)) +
  geom_text(aes(label=V1), vjust = -0.2)+
  ylab("N") +
  facet_wrap(~level, scales = "free_y")
ggsave("01_IDs_barplot.pdf", width = (nruns/2)+8)
dev.off()

# Proteus in-depth analysis
###########################

# Plot heatmaps of peptide and protein intensities as summarized by Proteus
evi = readEvidenceFile(file = evidenceFile)
pepdat <- makePeptideTable(evi, meta)

if(quantile_normalization == TRUE){
  # backup
  pepdat.raw = pepdat
  
  # Normalize using R/Limma quantile normalisation
  pepdat.n_quantile = normalizeData(pepdat, norm.fun = limma::normalizeQuantiles)
  protdat.n_quantile = makeProteinTable(pepdat.n_quantile)
  
  # pdf("03_Normalization.pdf", width = 3+nruns/3)
  boxplot(pepdat$tab, log="y", main = "Peptide intensities, Before normalization", las = 2)
  boxplot(pepdat.n_quantile$tab, log="y", main = "Peptide intensities, After quantile normalization, Limma", las = 2)
  # dev.off()
  
  pepdat = pepdat.n_quantile
  
}

pepmatrix = pepdat$tab
pepmatrix[is.na(pepmatrix)] = 0
pepmatrix = pepmatrix[rowSums(pepmatrix)>0,]
pepmatrix.log10 = log10(pepmatrix)
pepmatrix.log10[is.infinite(pepmatrix.log10)] = 0

# assemble annotation
annotation_column = pepdat$metadata
rownames(annotation_column) = annotation_column$experiment
annotation_column$experiment = NULL
annotation_column$sample = NULL
annotation_column$measure = NULL

pheatmap(pepmatrix.log10,
         # scale = "row",
         cluster_cols = T,
         show_rownames = F,
         ann_col = annotation_column,
         main = "Log10 Peptide intensities",
         color = viridis::inferno(30))

dev.copy(pdf, file = "03a_Heatmap_peptides.pdf",
         width = (nruns/2)+8,
         height = (nruns/2)+8)
dev.off()

# Plot heatmap of protein intensities as summarized by Proteus

# make Protein Table
protdat <- makeProteinTable(pepdat)
protdat$tab[is.na(protdat$tab)] = 0
protmatrix = protdat$tab[rowSums(protdat$tab)>0,]
protmatrix.log10 = log10(protmatrix)
protmatrix.log10[is.infinite(protmatrix.log10)] = 0

colnames(protmatrix.log10)
rownames(annotation_column)

pheatmap(protmatrix.log10,
         # scale = "row",
         cluster_cols = T,
         ann_col = annotation_column,
         show_rownames = F,
         ann_col = annotation_column,
         main = "Log10 Protein intensities (top3sum)",
         color = viridis::inferno(30))

dev.copy(pdf, file = "03b_Heatmap_proteins_top3sum.pdf",
         width = (nruns/2)+8,
         height = (nruns/2)+8)
dev.off()

# Plot Pearson correlation matrices
# Peptide level
plotDistanceMatrix(pepdat) + geom_text(aes(x = Sample1, y = Sample2, label = round(value,2))) +
  ggtitle("Peptide intensity correlation, Pearson") +
  scale_fill_gradient(low = "white", high = "red")
ggsave("05_Correlation_matrix_peptides.pdf", height = (nruns/3)+4, width = 4+(nruns/2.5))
dev.off()

# Protein level
plotDistanceMatrix(protdat) + geom_text(aes(x = Sample1, y = Sample2, label = round(value,2))) +
  ggtitle("Protein intensity correlation, Pearson") +
  scale_fill_gradient(low = "white", high = "red")
ggsave("05_Correlation_matrix_proteins.pdf", height = (nruns/3)+4, width = 4+(nruns/2.5))
dev.off()

# PairsPanel with smoothed density of xy scatter
for (condition in protdat$conditions){
  matrix = protmatrix.log10[, protdat$conditions %in% condition]
  matrix.completeobs = matrix[rowSums(matrix == 0) == 0,]
  pairs.panels(matrix.completeobs, smoother = TRUE, ellipses = F, lm = TRUE)
  dev.copy(pdf, paste0("06_correlation_pairpanel_", condition, ".pdf"), height = ncol(matrix)*4, width = ncol(matrix)*4)
  dev.off()
}

# Mass errors
ev2 = copy(ev)
names(ev2) = gsub("\\]|\\[|\\)|\\(|\\/| ", "_", names(ev2))
ggplot(ev2) + geom_density(aes(x = Uncalibrated_mass_error__ppm_, color = condition, group = Raw_file)) +
  ggtitle("Uncalibrated mass errors ppm") +
  theme_minimal()
ggsave("MassErrors.pdf", height = 5, width = 5)
