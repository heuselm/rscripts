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
evidenceFile = "../txt/evidence.txt"
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

# Plot IDs
nruns = length(unique(ev$`Raw file`))

ev = merge(ev, meta, by.x = "Experiment", by.y = "experiment")

ids.prot = ev[, length(unique(Proteins)), .(`Raw file`,Experiment,sample,condition,replicate)]
ids.pep = ev[, length(unique(`Modified sequence`)), .(`Raw file`,Experiment,sample,condition,replicate)]

ids = rbind(ids.prot[, level:="Protein.groups"],
            ids.pep[, level:="Modified.peptides"])

ggplot(ids, aes(x = Experiment, y = V1, fill = condition)) + geom_bar(stat = "identity") + ggtitle("Identifications MaxQuant") +
  theme(axis.text.x = element_text(angle = 90)) +
  geom_text(aes(label=V1), vjust = -0.2)+
  ylab("N") +
  facet_wrap(~level, scales = "free_y")
ggsave("01_IDs_manual.pdf", width = (nruns/2)+8)

ggplot(ev, aes(log10(Intensity), col = Experiment)) + geom_density() + facet_wrap(~condition, ncol = 1)
ggsave("01_IntensityDist_manual.pdf", height = (nruns/2)+8)

# Proteus differential analysis
################################
evi = readEvidenceFile(evidenceFile)

pepdat <- makePeptideTable(evi, meta)
protdat <- makeProteinTable(pepdat)

# Plot ID numbers
plotCount(pepdat)
ggsave("02_pepCount_proteus.pdf")

plotCount(protdat)
ggsave("02_protCount_proteus_min2pep.pdf")

# Normalize using R/Limma quantile normalisation
pepdat.n_quantile = normalizeData(pepdat, norm.fun = limma::normalizeQuantiles)
protdat.n_quantile = makeProteinTable(pepdat.n_quantile)

pdf("03_Normalization.pdf", width = 3+nruns/3)
boxplot(pepdat$tab, log="y", main = "Peptide intensities, Before normalization", las = 2)
boxplot(pepdat.n_quantile$tab, log="y", main = "Peptide intensities, After quantile normalization, Limma", las = 2)
dev.off()

pdf("03_Normalization_pep_impact_on_prot.pdf", width = 3+nruns/3)
boxplot(protdat$tab, log="y", main = "Protein intensities, Before pep-level normalization", las = 2)
boxplot(protdat.n_quantile$tab, log="y", main = "Protein intensities, After pep-level quantile normalization, Limma", las = 2)
dev.off()

# Count common vs total ids, total dist. of all pairwise comp. (Jaccard Index)
plotDetectionSimilarity(pepdat.n_quantile, bin.size = 0.02) + ggtitle("Jaccard Similarity")
ggsave("04_detectionSimilarity.pdf")

# Plot Pearson correlation matrices
# Peptide level
plotDistanceMatrix(pepdat) + geom_text(aes(x = Sample1, y = Sample2, label = round(value,2))) +
  ggtitle("Peptide intensity correlation, Pearson, raw peptide intensities") +
  scale_fill_gradient(low = "white", high = "red")
ggsave("05_correlation_matrix_peptides_raw.pdf", height = nruns/3, width = nruns/2.5)
plotDistanceMatrix(pepdat.n_quantile) + geom_text(aes(x = Sample1, y = Sample2, label = round(value,2))) +
  ggtitle("Peptide intensity correlation, Pearson, normalized peptide intensities") +
  scale_fill_gradient(low = "white", high = "red")
ggsave("05_correlation_matrix_peptides_norm.pdf", height = nruns/3, width = nruns/2.5)

# Protein level
plotDistanceMatrix(protdat) + geom_text(aes(x = Sample1, y = Sample2, label = round(value,2))) +
  ggtitle("Protein intensity correlation, Pearson, using raw peptide intensities") +
  scale_fill_gradient(low = "white", high = "red")
ggsave("05_correlation_matrix_proteins_raw.pdf", height = nruns/3, width = nruns/2.5)
plotDistanceMatrix(protdat.n_quantile) + geom_text(aes(x = Sample1, y = Sample2, label = round(value,2))) +
  ggtitle("Protein intensity correlation, Pearson, using normalized peptide intensities") +
  scale_fill_gradient(low = "white", high = "red")
ggsave("05_correlation_matrix_proteins_norm.pdf", height = nruns/3, width = nruns/2.5)

# Plot some example peptides/proteins
pdf("06_protein_plots_rand10test.pdf", width = nruns/3)
for (i in 1:20){
  pi = plotProtPeptides(pepdat = pepdat.n_quantile, protein = sample(pepdat$proteins, 1))
  plot(pi)
}
dev.off()

# Principal component analysis
pdf("07_PCA_plots_separate.pdf", height = 6, width = 7)
plotPCA(pepdat, point.size = 3, label.size = 1) + ggtitle("Peptide data, raw")
plotPCA(pepdat.n_quantile, point.size = 3, label.size = 1) + ggtitle("Peptide data, quantile-normalised")
plotPCA(protdat, point.size = 3, label.size = 1) + ggtitle("Protein data, raw")
plotPCA(protdat.n_quantile, point.size = 3, label.size = 1) + ggtitle("Protein data, quantile-normalised")
dev.off()

pca1 = plotPCA(pepdat, point.size = 3, label.size = 1) + ggtitle("Peptide data, raw")
pca2 = plotPCA(pepdat.n_quantile, point.size = 3, label.size = 1) + ggtitle("Peptide data, quantile-normalised")
pca3 = plotPCA(protdat, point.size = 3, label.size = 1) + ggtitle("Protein data, raw")
pca4 = plotPCA(protdat.n_quantile, point.size = 3, label.size = 1) + ggtitle("Protein data, quantile-normalised")

pcas = gridExtra::grid.arrange(pca1,pca2,pca3,pca4, top = "Principal Component Analysis")
ggsave("07_PCA_plots_combined.pdf", plot = pcas, height = 12, width = 14)

# clustering
clust_1 = plotClustering(pepdat) + ggtitle("based on peptide data, raw")
clust_2 = plotClustering(pepdat.n_quantile) + ggtitle("based on peptide data, normalized")
clust_3 = plotClustering(protdat) + ggtitle("based on protein data, raw")
clust_4 = plotClustering(protdat.n_quantile) + ggtitle("based on protein data, normalized")

Dendrograms = gridExtra::grid.arrange(clust_1,clust_2,clust_3,clust_4,
                        top = "Sample-Sample clustering Analysis (Pearson correlation, complete obs.)")
ggsave("08_Clustering_Dendrograms_Pearson.pdf", Dendrograms, height = 12, width = 12)

# Protein level Mean-variance relationship
mv1 = plotMV(pepdat, with.loess = T) + ggtitle("Peptide data, raw")
mv2 = plotMV(pepdat.n_quantile, with.loess = T) + ggtitle("Peptide data, quantile-normalised")
mv3 = plotMV(protdat, with.loess = T) + ggtitle("Protein data, raw")
mv4 = plotMV(protdat.n_quantile, with.loess = T) + ggtitle("Protein data, quantile-normalised")

MVs = gridExtra::grid.arrange(mv1,mv2,mv3,mv4,
                        top = "Variance over mean log10 intensity relationship analysis")
ggsave("09_Mean_Variance_Relationships.pdf", MVs, height = 14, width = 14)


## Differential expresssion testing between (2) conditions
# all possible pairwise comparisons
pairs = t(combn(unique(meta$condition), 2))
pairs

plot_adjusted_pvals = TRUE
plot_protpep_plots = FALSE

# # impute missing values
# library(pheatmap)
# tab = 
# # View(tab)
tab_log10 = log10(protdat.n_quantile$tab)
# tab_log10[is.na(tab_log10)] = 0

tab_log10[is.na(tab_log10)] = rnorm(n = length(tab_log10[!is.na(tab_log10)]), mean = quantile(tab_log10[!is.na(tab_log10)], 0.005), sd = 0.2)
# View(tab_log10)
tab_imp = 10^tab_log10
protdat.n_quantile.imp = copy(protdat.n_quantile)
protdat.n_quantile.imp$tab = tab_imp

library(pheatmap)
pheatmap(tab_log10, na.rm = TRUE)
dev.copy(pdf, file = "Heatmap_log10int_all_imp.pdf", height = 15, width = 7)
dev.off()

for (i in 1:nrow(pairs)){
  c1 = pairs[i,1]
  c2 = pairs[i,2]
  
  diri = paste0("diffTest_",c2,"vs",c1, "_imp")
  dir.create(diri)
  
  res <- limmaDE(protdat.n_quantile.imp, formula = "~condition", conditions = pairs[i,], 
               transform.fun = log2, sig.level = 0.05)
  names(res)
  if(plot_adjusted_pvals){
    res$P.Value.backup = res$P.Value
    res$P.Value = res$adj.P.Val
  }
  
  p1 = plotVolcano(res, binhex = F, plot.grid = T) +
  ggtitle(paste("Left: Higher abundant in", c1, "   ", "Right: Higher abundant in", c2)) +
  geom_hline(yintercept = -log10(0.05), col = "red", lty = 2) +
    geom_vline(xintercept = c(-log2(2),log2(2)), , col = "red", lty = 2) +
    # scale_x_continuous(breaks = seq(-10,10,1)) +
    geom_text_repel(aes(label = protein))
  p1
  ggsave(paste0(diri,"/differential_expression_Volcanoplot.pdf"), height = 8, width = 8)
  
  # p2 = plotVolcano(res, binhex = F, marginal.histograms = T, plot.grid = T)
  # p2
  # ggsave(paste0(diri,"/differential_expression_Volcanoplot_02.pdf", p2))
  
  htmlwidgets::saveWidget(ggplotly(p1 + geom_point(aes(text = protein, alpha = -log10(P.Value)))),
                          file = paste0(getwd(),"/",diri,"/differential_expression_Volcanoplot_interactive.html"))
  
  res$comparison = paste0(c2,"vs",c1)
  fwrite(res, file = paste0(diri,"/differential_expression_table.tsv"), sep = "\t")
  
  ## Plot significant proteins intensities
  #c1
  pdf(paste0(diri,"/differential_proteins_higher_in_",c1,".pdf"), width = nruns/3)
  
  res_prots_up_c1 = na.omit(res[res$significant == TRUE & res$logFC < 0,])
  res_prots_up_c1 = res_prots_up_c1[order(res_prots_up_c1$logFC),]
  
  # include highlighted volcano
  plot(p1 + geom_point(data = res_prots_up_c1, size = 3, col = "red") +
    geom_text_repel(data = res_prots_up_c1, aes(label = protein), force = 10, max.iter = 50000, max.overlaps = 10000) +
    ggtitle(paste("Proteins higher abundant in", c1)))
  
  if(plot_protpep_plots){
    prots_up_c1 = res_prots_up_c1$protein
    for (i in seq_along(prots_up_c1)){
      pi = plotProtPeptides(pepdat = pepdat.n_quantile, protein = prots_up_c1[i])
      plot(pi)
    }  
  }
  
  
  dev.off()
  
  # c2
  pdf(paste0(diri,"/differential_proteins_higher_in_",c2,".pdf"), width = nruns/3)
  res_prots_up_c2 = na.omit(res[res$significant == TRUE & res$logFC > 0,])
  res_prots_up_c2 = res_prots_up_c2[order(-res_prots_up_c2$logFC),]
  
  # include highlighted volcano
  plot(p1 + geom_point(data = res_prots_up_c2, size = 3, col = "red") +
    geom_text_repel(data = res_prots_up_c2, aes(label = protein), force = 10, max.iter = 50000, max.overlaps = 10000) +
    ggtitle(paste("Proteins higher abundant in", c2)))
  
  prots_up_c2 = res_prots_up_c2$protein
  if(plot_protpep_plots){
  for (i in seq_along(prots_up_c2)){
    pi = plotProtPeptides(pepdat = pepdat.n_quantile, protein = prots_up_c2[i])
    plot(pi)
    }
  }
  dev.off()
  
}

res_M1vsGFP = fread("diffTest_M1_0mMvsGFP_0mM_imp/differential_expression_table.tsv")

res_M1vsGFP[, M1_binder:=FALSE]
res_M1vsGFP[P.Value <= 0.01 & logFC > 2, M1_binder:=TRUE]

quant_significant_proteins = tab_log10[ row.names(tab_log10) %in% res_M1vsGFP[P.Value <= 0.01 & logFC > 2]$protein, ]

pheatmap(tab_log10[ row.names(tab_log10) %in% res_M1vsGFP[P.Value <= 0.01 & logFC > 2]$protein, ])
dev.copy(pdf, file = "Heatmap_log10int_M1binders_imp.pdf", height = 7, width = 7)
dev.off()



