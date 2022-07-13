# most simple test of running a plotting script with Rscript and user selection of input file without touching R
if(is.null(require('mocode'))){
devtools::install_github("heuselm/mocode")
}

# Select inputfile/folder
inputfile <- choose.files(caption = "Place DIANN results (ResA.tsv, ResB.tsv and ResA.stats.tsv, ResB.stats.tsv) into one folder and choose one of these files. All resultsets in the folder will be compared.")

# plot comparative plots
mocode::compareDiannResultSets(result_folder = dirname(inputfile),
Venn_diagram = TRUE,
compare_precursor_level = FALSE,
heatmap = TRUE,
sample_correlation = TRUE)
