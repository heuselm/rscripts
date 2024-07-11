#' ---
#' title: "Summary Report of DIA-NN resultset"
#' author: mhe
#' output:
#'   html_document:
#'     keep_tex: false
#' ---

#'## What is this?
#'Summary of DIANN output table.tsv (long format) with consideration of the study design.

#'## Environment setup

#+ r load packages
library('data.table')
library('ggplot2')
library('pheatmap')
library('plotly')
library('corrplot')
library('psych')

# custom color scheme

#'## Data import & annotation with study design (condition & replicate)
#+ Data import
# path_to_dataset = "S:\\Applications\\P003_500SPD\\03_Results\\220705_ShortCycleDiaMethodTest\\report_LibFree_2.tsv"
if (!exists("path_to_dataset")) {
  path_to_dataset <- choose.files(caption = "Choose report.tsv file from DIANN output folder")
  data = fread(path_to_dataset)
} 

setwd(dirname(path_to_dataset))
dsname = basename(path_to_dataset)

# write template file for user
nruns = length(unique(data$File.Name))
fwrite(data.table("File.Name" = unique(data$File.Name),
                  "Condition" = rep("FILL", nruns),
                  "Replicate" = rep("FILL", nruns)),
       file = "StudyDesign_template.csv")

message('StudyDesign_template.csv written to result folder - fill in condition and replicate info and save as StudyDesign.csv')

# read & merge StudyDesign
study_design = fread('StudyDesign.csv')

data = merge(data, study_design, by = 'File.Name', all.x = F)
data[, Replicate:=factor(Replicate, levels = sort(unique(Replicate)))]
data[, Condition:=factor(Condition, levels = sort(unique(Condition)))]

# IDs, per run
#+ IDcount per run
ggplot(data[, length(unique(Precursor.Id)), .(Condition,Replicate,File.Name)],
       aes(Condition, group = File.Name, y = V1)) +
    geom_bar(stat = "identity", position = position_dodge2(preserve = 'single'),
             width = 0.8, col = "black") +
     ylab("N Precursors") +
  # scale_fill_manual(values = customcolors[c(3,1)]) +
    ggtitle("Precursors Identified") +
  theme_minimal() +
   theme(axis.text.x = element_text(angle = 45, hjust = 1)) 
ggsave(paste0("IDs_Precursors_",dsname,".pdf"), width = 5+(nruns/6), height = 6)

ggplot(data[, length(unique(Protein.Group)), .(Condition,Replicate,File.Name)],
       aes(Condition, group = File.Name, y = V1)) +
    geom_bar(stat = "identity", position = position_dodge2(preserve = 'single'),
             width = 0.8, col = "black") +
     ylab("N Protein Groups") +
  # scale_fill_manual(values = customcolors[c(3,1)]) +
    ggtitle("Protein Groups Identified") +
  theme_minimal() +
   theme(axis.text.x = element_text(angle = 45, hjust = 1)) 
ggsave(paste0("IDs_ProteinGroups_",dsname,".pdf"), width = 12, height = 6)


# IDs, mean and average
#+ IDcount mean and average, fig.width=12
# Now, we calculate average IDs and protein intensity CVs
IDs_prec = data[, length(unique(Precursor.Id)), .(Condition, Replicate)]
IDs_prot = data[, length(unique(Protein.Group)), .(Condition, Replicate)]
IDs = rbind(IDs_prec[, level:="Precursors"],
            IDs_prot[, level:="Protein.Groups"])
setnames(IDs, "V1", "N")

# Calc Mean + SD
IDs[, mean:=mean(N), .(Condition, level)]
IDs[, sd:=sd(N), .(Condition, level)]
IDs[, N:=NULL]
IDs[, Replicate:=NULL]
IDs = unique(IDs)


# Plot
ggplot(IDs[level == "Protein.Groups"], aes(x = Condition, y = mean)) +
 geom_bar(stat = "identity") +
 geom_errorbar(aes(ymin = mean-sd, ymax = mean+sd), width = 0.2) +
  theme_minimal() +
 theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
 # scale_fill_manual(values = customcolors[1]) +
 ggtitle("Protein.Groups identified") 
ggsave(paste0("IDs_ProteinGroups_mean_sd_",dsname,".pdf"), width = 5+(nruns/12), height = 5)

ggplot(IDs[level == "Precursors"], aes(x = Condition, y = mean)) +
 geom_bar(stat = "identity") +
 geom_errorbar(aes(ymin = mean-sd, ymax = mean+sd), width = 0.2) +
 theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme_minimal() +
 # scale_fill_manual(values = customcolors[1]) +
 ggtitle("Precursors identified")
ggsave(paste0("IDs_Precursors_mean_sd_", dsname, ".pdf"), width = 5+(nruns/12), height = 5)

#'## Protein quantification Coefficients of Variation (CV analysis, PG.MaxLFQ)
#+ ProteinQuantCV

# define coefficient of variation function
cv = function(x){sd(x)/mean(x)*100}

# Calculate number of observations per condition
data[, PG.MaxLFQ.nobs:=0]
data[PG.MaxLFQ > 0, PG.MaxLFQ.nobs:=length(unique(Replicate)), .(Protein.Group, Condition)]

# calculate CV for those with n >= 2 and plot CV dists, PG.Normalised
data[PG.MaxLFQ.nobs >= 2, PG.MaxLFQ.CV:=cv(PG.MaxLFQ), .(Condition, Protein.Group)]

data_protCv = unique(data[, .(Condition,Protein.Group,PG.MaxLFQ.CV)])


# Density plot, count-scaled y axis
ggplot(data_protCv, aes(x = PG.MaxLFQ.CV, y = ..count.., col = Condition)) +
  geom_density(size = 1) +
  scale_color_brewer(palette = "Set1") + 
  theme_minimal() +
  coord_cartesian(xlim = c(0,50)) +
  ggtitle("Protein quantification CV, density")
ggsave(paste0("CVs_PG.MaxLFQ_n3_density_nscaled_", dsname,".pdf"), height = 7, width = 7)

# Density plot, un-scaled y axis
ggplot(data_protCv, aes(x = PG.MaxLFQ.CV, col = Condition)) +
  geom_density(size = 1) +
  scale_color_brewer(palette = "Set1") + 
  theme_minimal() +
  coord_cartesian(xlim = c(0,50)) +
  ggtitle("Protein quantification CV, density")
ggsave(paste0("CVs_PG.MaxLFQ_n3_density_", dsname, ".pdf"), height = 7, width = 7)

# Violin & Box-plot
ggplot(data_protCv, aes(x = Condition, y = PG.MaxLFQ.CV)) + 
  geom_violin(aes(fill = Condition)) +
  geom_boxplot(width = 0.1, outlier.shape = NA) +
  scale_color_brewer(palette = "Set1") + 
  scale_fill_brewer(palette = "Set1") + 
  theme_minimal() +
  coord_cartesian(ylim = c(0,50)) +
  ggtitle("Protein MaxLFQ quantification CVs") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(paste0("CVs_PG.MaxLFQ_n3_violin-box_", dsname, ".pdf"), height = 7, width = 14)

#'## Intensity distributions
#+ Intensity distributions
# Precursor level intensities, RAW intensities
ggplot(data) + geom_violin(aes(x = Condition,
                                y = log10(Precursor.Quantity),
                                fill = Condition, lty = Replicate), alpha = 0.2,
                            position = position_dodge(width=0.8)) +
  geom_boxplot(aes(x = Condition,
                  y = log10(Precursor.Quantity),
                  lty = Replicate), width = 0.5,
              position = position_dodge(width=0.8)) +
  scale_fill_brewer(palette = "Set1") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle("Precursor.Quantity distribution RAW")
ggsave(paste0("IntensityDistribution_Precursor.Quantity_Violin_", dsname, ".pdf"), height = 5, width = 5+(nruns/12))

# Precursor level intensities, Normalised intensities
ggplot(data) + geom_violin(aes(x = Condition,
                               y = log10(Precursor.Normalised),
                               fill = Condition, lty = Replicate), alpha = 0.2,
                           position = position_dodge(width=0.8)) +
  geom_boxplot(aes(x = Condition,
                   y = log10(Precursor.Normalised),
                   lty = Replicate), width = 0.5,
               position = position_dodge(width=0.8)) +
  scale_fill_brewer(palette = "Set1") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle("Precursor.Normalised distribution")
ggsave(paste0("IntensityDistribution_Precursor.Normalised_Violin_", dsname, ".pdf"), height = 5, width = 5+(nruns/12))


#'## IDs per Minute Gradient time
#+IDsPerMinute
# For all datasets:
data[, RT_bin_minutes:=cut(RT, seq(0,round(max(data$RT), 1)))]
idpm = unique(data[, length(unique(Precursor.Id)),.(Condition,Replicate,File.Name, RT_bin_minutes)])

ggplot(idpm) + geom_bar(aes(x = RT_bin_minutes, y = V1, fill = Condition),
                             stat = "identity", position = position_dodge2(preserve = "single")) +
  xlab("Retention time bin") + ylab("Precursor identifications per minute") +
   ggtitle("Precursor identifications per Retention time bin, minute") + scale_fill_brewer(palette = "Set1") +
  theme_minimal()
ggsave(paste0("IDs_PrecursorsPerMinute_bar_", dsname, ".pdf"), height = 5, width = 5+(nruns/12))


#'## Pairwise log-log correlation analysis
#+ Pairwise log-log correlation analysis
data[, Condition:=factor(Condition)]
data[, Replicate:=as.numeric(Replicate)]

pdf("ReplicateCorrelation_log10MaxLFQ_R1R2.pdf")
for (condition in unique(data$Condition)){
      m.dt = dcast(unique(data[Condition == condition & Replicate <= 2,
                               .(Protein.Group, PG.MaxLFQ, Condition, Replicate)]),
                   Protein.Group~Condition+Replicate, value.var = "PG.MaxLFQ")
      m = as.matrix(m.dt[,2:3])
      rownames(m) = m.dt$Protein.Group
      colnames(m) = paste(condition, c("R1","R2"))
      m.log10 = log10(m)
      m.log10[is.na(m.log10)] = 0
      m.log10[is.infinite(m.log10)] = 0
      m.nobs = 2-rowSums(is.na(m))
     
      pairs.panels(m.log10[m.nobs == 2,], smoother = TRUE, lm = TRUE)
}
dev.off()

pdf(paste0("ReplicateCorrelation_log10MaxLFQ_all_", dsname, ".pdf"))
for (condition in unique(data$Condition)){
  m.dt = dcast(unique(data[Condition == condition,
                           .(Protein.Group, PG.MaxLFQ, Condition, Replicate)]),
               Protein.Group~Condition+Replicate, value.var = "PG.MaxLFQ")
  m = as.matrix(m.dt[,2:ncol(m.dt)])
  rownames(m) = m.dt$Protein.Group
  colnames(m) = paste(condition, paste0("R", seq(1, ncol(m))))
  m.log10 = log10(m)
  m.log10[is.na(m.log10)] = 0
  m.log10[is.infinite(m.log10)] = 0
  m.nobs = ncol(m)-rowSums(is.na(m))
  
  pairs.panels(m.log10[m.nobs == ncol(m),], smoother = TRUE, lm = TRUE)
}
dev.off()

# ID density / 'IDTic'
# names(data)
ggplot(data[!grepl("blank", Condition)]) +
  geom_density(aes(x = RT, y = ..count.., color = Condition, group = File.Name), adjust = 0.2) +
  theme_minimal() + 
  ggtitle("Identification density by Condition")
ggsave(paste0("Identification density by Condition_", dsname, ".pdf"), height = 5, width = 10)

ggplot(data[!grepl("blank", Condition)]) +
  geom_density(aes(x = RT, y = ..count.., color = Condition, group = File.Name), adjust = 0.2) +
  theme_minimal() + facet_wrap(~Condition) +
  ggtitle("Identification density by Condition")
ggsave(paste0("Identification density by Condition split by Condition", dsname, ".pdf"), height = 5, width = 15)
