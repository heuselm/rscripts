#' Plot MaxQuant QC Summary
#' QC of yeast injections on Thermo OT HFX 'Minerva'
#####################################################
# load packages
library(data.table)
library(ggplot2)

#' First, plot summary metrics
# load n prep data
d = fread("summary.txt")
d = d[1:5]
setorder(d, Experiment)
names(d)

# plot
barplot(d$MS, beside = TRUE, names.arg = d$Experiment,
        main = "MS scans")
barplot(d$Peaks, beside = TRUE, names.arg = d$Experiment,
        main = "Peaks")
barplot(d$`Peaks Sequenced [%]`, beside = TRUE, names.arg = d$Experiment,
        main = "Peaks Sequenced [%]")
barplot(d$`MS/MS`, beside = TRUE, names.arg = d$Experiment,
        main = "MS/MS scans")

barplot(d$`MS/MS Submitted`, beside = TRUE, names.arg = d$Experiment,
        main = "MSMS submitted")
barplot(d$`MS/MS Identified`, beside = TRUE, names.arg = d$Experiment,
        main = "MSMS identified")
barplot(d$`Peptide Sequences Identified`, beside = TRUE, names.arg = d$Experiment,
        main = "Peptide Sequences")
barplot(d$`Av. Absolute Mass Deviation [ppm]`, beside = TRUE, names.arg = d$Experiment,
        main = "`Av. Absolute Mass Deviation [ppm]`")



# output as separate .pdf
pdf("summary.pdf")
par(mfrow = c(4,2))
barplot(d$MS, beside = TRUE, names.arg = d$Experiment,
        main = "MS scans")
barplot(d$Peaks, beside = TRUE, names.arg = d$Experiment,
        main = "Peaks")
barplot(d$`Peaks Sequenced [%]`, beside = TRUE, names.arg = d$Experiment,
        main = "Peaks Sequenced [%]")
barplot(d$`MS/MS`, beside = TRUE, names.arg = d$Experiment,
        main = "MS/MS scans")

barplot(d$`MS/MS Submitted`, beside = TRUE, names.arg = d$Experiment,
        main = "MSMS submitted")
barplot(d$`MS/MS Identified`, beside = TRUE, names.arg = d$Experiment,
        main = "MSMS identified")
barplot(d$`Peptide Sequences Identified`, beside = TRUE, names.arg = d$Experiment,
        main = "Peptide Sequences")
barplot(d$`Av. Absolute Mass Deviation [ppm]`, beside = TRUE, names.arg = d$Experiment,
        main = "`Av. Absolute Mass Deviation [ppm]`")
dev.off()

#' Second, summarize MSMS information to try and shed light on why less msms led to IDs
msmsScans = fread("msmsScans.txt")
names(msmsScans) = gsub(" ", "_", names(msmsScans))
setorder(msmsScans, "Experiment")
msmsScans[, Raw_file:=factor(Raw_file, levels = unique(Raw_file))]

# MSMS TIC violin
ggplot(msmsScans, aes(Raw_file, as.numeric(Total_ion_current))) +
  geom_violin() +
  theme_classic() +
  geom_violin(aes(y=as.numeric(Base_peak_intensity), fill = "red", alpha = 0.3)) +
  scale_y_log10() + theme(legend.position = "none") +
  ggtitle("MSMS TIC and BP (red) intensities")

# MSMS TIC
ggplot(msmsScans, aes(Retention_time, as.numeric(Total_ion_current))) +
  geom_line() + scale_color_brewer(palette = "Spectral") +
  theme_classic() + facet_wrap(~Raw_file, nrow = 1) +
  ggtitle("MSMS TIC")

# MSMS Fill times
ggplot(msmsScans, aes(Retention_time, Ion_injection_time)) +
  geom_point(pch = ".") + scale_color_brewer(palette = "Spectral") +
  theme_classic() + facet_wrap(~Raw_file, nrow = 1) +
  ggtitle("MSMS fill times")

# MSMS identifications over RTnMZ
ggplot(msmsScans, aes(Retention_time, `m/z`, alpha = Identified)) +
  geom_point(pch = ".")  +
  theme_classic() + facet_wrap(~Raw_file, nrow = 1) +
  ggtitle("MSMS identified mzrt")

#' MS2 TIC dropped while fill times remain stable, i.e. AGC control is confused / ion transmission for MSMS seems
#' to be impaired.


#' Look at mass error distribution from evidence.txt table
evidence = fread("evidence.txt")
names(evidence) = gsub(" ", "_", names(evidence))
names(evidence) = gsub("\\[", "_", names(evidence))
names(evidence) = gsub("\\]", "_", names(evidence))
names(evidence)
setorder(evidence, Experiment)
evidence[, Raw_file:=factor(Raw_file, levels = unique(Raw_file))]

# Mass errors
ggplot(evidence, aes(x = Raw_file, y = Uncalibrated_Mass_Error__ppm_)) + geom_violin() +
  geom_violin(aes(y = Mass_Error__ppm_), fill = "blue", alpha = 0.3) +
  ylab("Mass Error ppm uncal/cal") + ggtitle("Mass Errors uncalibrated vs calibrated")

#' Mass accuracy is drifting a bit but not too much.
#' Service required.