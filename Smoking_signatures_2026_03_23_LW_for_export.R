# Smoking mutational signatures
# 2023-12-10 - updated by LWoodhouse 2026-04-08 for export

#Load libraries and set wd ------------------------------------------------------
#Smoking mutational signatures
setwd("~/personal_file_path/")
options(stringsAsFactors=F)
library(tidyverse)
library(dplyr)
library(reshape2)
library(ggsignif)
library(cowplot)
library(ggpubr)
library(eulerr)
library(survival)
library(survminer)
library(VennDiagram)
library(gridExtra)
library(rstatix)

par(mar=c(1,4,1,1))

# Import and process dataset -------------------------------------------------------
#Cohort using the v16 release and the following filters: 
#Disease type equals "LUNG", 
#Disease subtype equals one of "ADENOCARCINOMA; SQUAMOUS_CELL; ADENOSQUAMOUS"
#Preparation method equals "FF"
#Somatic coding variants per Mb > 0.03 (excludes samples with <100 SNVs)
#Study abbreviation equals one of: LUG; LUAD; LUSC
#Produces 980 NSCLC participants

#Import and process clinical data
clinical_data <- read.csv("~/file/path/100K_lung_cancer_clinical_data_updated.csv", header = T)
clinical_data[clinical_data==""]  <- NA 
clinical_data <- clinical_data[,-c(25:60)]

#Import and process signatures 
#SBS signatures
#Process and filter data
COSMIC_SBS_all <- read.table(file = "/file/path/Activities/COSMIC_SBS96_Activities_refit.txt", header = T)
COSMIC_SBS_all$Samples <- sub("[_].*", "", COSMIC_SBS_all$Samples)
COSMIC_SBS <- COSMIC_SBS_all[COSMIC_SBS_all$Samples %in% clinical_data$Gel_ID, ]
COSMIC_SBS <- COSMIC_SBS[, colSums(COSMIC_SBS != 0) !=0]
rownames(COSMIC_SBS) <- COSMIC_SBS$Samples
COSMIC_SBS <- COSMIC_SBS[,-1]
# Calculate percentages and clean up
COSMIC_SBS_percentages <- COSMIC_SBS/rowSums(COSMIC_SBS)
COSMIC_SBS_percentages[COSMIC_SBS_percentages == 0] <- NA
SBS <- COSMIC_SBS_percentages[, colSums(!is.na(COSMIC_SBS_percentages))>0]

#DBS signatures
#Process and filter data
COSMIC_DBS_all <- read.table(file = "/file/path/Activities/COSMIC_DBS78_Activities_refit.txt", header = T)
COSMIC_DBS_all$Samples <- sub("[_].*", "", COSMIC_DBS_all$Samples)
COSMIC_DBS <- COSMIC_DBS_all[COSMIC_DBS_all$Samples %in% clinical_data$Gel_ID, ]
COSMIC_DBS <- COSMIC_DBS[, colSums(COSMIC_DBS != 0) !=0]
rownames(COSMIC_DBS) <- COSMIC_DBS$Samples
COSMIC_DBS <- COSMIC_DBS[,-1]
#Calculate percentages and clean up
COSMIC_DBS_percentages <- COSMIC_DBS/rowSums(COSMIC_DBS)
COSMIC_DBS_percentages[COSMIC_DBS_percentages == 0] <- NA
DBS <- COSMIC_DBS_percentages[, colSums(!is.na(COSMIC_DBS_percentages))>0]

#ID signatures
#Process and filter data 
COSMIC_IDS_all <- read.table(file = "f/ile/path/Activities/COSMIC_ID83_Activities_refit.txt", header = T)
COSMIC_IDS_all$Samples <- sub("[_].*", "", COSMIC_IDS_all$Samples)
COSMIC_IDS <- COSMIC_IDS_all[COSMIC_IDS_all$Samples %in% clinical_data$Gel_ID,]
COSMIC_IDS <- COSMIC_IDS[, colSums(COSMIC_IDS != 0) !=0]
rownames(COSMIC_IDS) <- COSMIC_IDS$Samples
COSMIC_IDS <- COSMIC_IDS[,-1]
#Calculate percentages and clean up
COSMIC_IDS_percentages <- COSMIC_IDS/rowSums(COSMIC_IDS)
COSMIC_IDS_percentages[COSMIC_IDS_percentages == 0] <- NA
IDS <- COSMIC_IDS_percentages[, colSums(!is.na(COSMIC_IDS_percentages))>0]

#Combine to create final signatures matrix with clinical details
SBS <- SBS %>% as.data.frame() %>% rownames_to_column("Gel_ID")
DBS <- DBS %>% as.data.frame() %>% rownames_to_column("Gel_ID")
IDS <- IDS %>% as.data.frame() %>% rownames_to_column("Gel_ID")
signatures <- purrr::reduce(.x = list(SBS, DBS, IDS), merge, by = "Gel_ID", all = TRUE) %>% replace(is.na(.), 0) 

lung_patients_GEL <- clinical_data[which(clinical_data$Gel_ID %in% signatures$Gel_ID),]
lung_patients_GEL <- purrr::reduce(.x = list(lung_patients_GEL, signatures), merge, by = "Gel_ID", all = TRUE)
lung_patients <- lung_patients_GEL

#Pseudoanonymised data with percentages
lung_patients_no_ID <- lung_patients_GEL
lung_patients_no_ID$Gel_ID <- paste("Patient_", 1:nrow(lung_patients_no_ID), sep = "")
SBS$Gel_ID <- lung_patients_no_ID$Gel_ID
DBS$Gel_ID <- lung_patients_no_ID$Gel_ID
IDS$Gel_ID <- lung_patients_no_ID$Gel_ID

#Pseudoanonymised signatures with raw counts 
SBS_raw_counts <- COSMIC_SBS
rownames(SBS_raw_counts) <- lung_patients_no_ID$Gel_ID
DBS_raw_counts <- COSMIC_DBS
rownames(DBS_raw_counts) <- lung_patients_no_ID$Gel_ID
IDS_raw_counts <- COSMIC_IDS
rownames(IDS_raw_counts) <- lung_patients_no_ID$Gel_ID

#Save signatures as CSV to export
write.csv(SBS_raw_counts, file = "~/file/path/Figures_final/SBS_raw_counts.csv")
write.csv(DBS_raw_counts, file = "~/file/path/Figures_final/DBS_raw_counts.csv")
write.csv(IDS_raw_counts, file = "~/file/path/Figures_final/IDS_raw_counts.csv")

#Pseudoanonymised data conversion key
conversion <- data.frame(cbind(lung_patients_GEL$Gel_ID,lung_patients_no_ID$Gel_ID))
write.csv(conversion, file = "~/file/path/Figures_final/Gel_ID_to_annonymised_key.csv")

#Combined clinical data with Gel ID and patient ID
lung_patients_combined <- data.frame(cbind(lung_patients_GEL,lung_patients_no_ID$Gel_ID))
write.csv(lung_patients_combined, file = "~/file/path/100k_lung_signatures_clinical_full.csv")

# Create cohorts for signatures
#SBS4 cohorts
SBS4 <- lung_patients_combined %>% filter(SBS4>0) #107
no_SBS4 <- lung_patients_combined %>% filter(SBS4==0) #22

#DBS2 cohorts
DBS2 <- lung_patients_combined %>% filter(DBS2>0) #122
no_DBS2 <- lung_patients_combined %>% filter(DBS2==0) #7

#ID3 cohorts
ID3 <- lung_patients_combined %>% filter(ID3>0) #109
no_ID3 <- lung_patients_combined %>% filter(ID3==0) #20

#Smoking cohorts
current_smoker <- lung_patients_combined %>% filter(Smoking.status == "Current") #31
ex_smoker <- lung_patients_combined %>% filter(Smoking.status == "Ex") #88
never_smoker <- lung_patients_combined %>% filter(Smoking.status == "Never") #8
ever_smoker <- lung_patients_combined %>% filter(Smoking.status %in% c("Current","Ex")) #119
missing_smoking_status <- lung_patients_combined %>% filter(Smoking.status = NA)

# SBS4 with smoking cohort
# Ever smokers and SBS4
ever_smoker_SBS4 <- ever_smoker %>% filter(SBS4>0) #104
ever_smoker_no_SBS4 <- ever_smoker %>% filter(SBS4 == 0) #15

# Ex-smokers and SBS4
ex_smoker_SBS4 <- ex_smoker %>% filter(SBS4>0) #73
ex_smoker_no_SBS4 <- ex_smoker %>% filter(SBS4 == 0) #15

# Never smokers and SBS4
never_smoker_SBS4 <- never_smoker %>% filter(SBS4>0) #1
never_smoker_no_SBS4 <- never_smoker %>% filter(SBS4 == 0) #7

# Sensitivity and specificity calculations for SBS4
sensitivity_SBS4 <- nrow(ever_smoker_SBS4)/nrow(ever_smoker) #0.87395
specificity_SBS4 <- nrow(never_smoker_no_SBS4)/nrow(never_smoker) #0.875
PLR_SBS4 <- sensitivity_SBS4/(1-specificity_SBS4) #6.991597
NLR_SBS4 <- (1-sensitivity_SBS4)/specificity_SBS4 #0.1440576
PPV_SBS4 <- nrow(ever_smoker_SBS4)/nrow(SBS4) #0.9719626
NPV_SBS4 <- nrow(never_smoker_no_SBS4)/nrow(no_SBS4) #0.3181818

# DBS2 with smoking cohort
# Ever smokers and DBS2
ever_smoker_DBS2 <- ever_smoker %>% filter(DBS2>0) #114
ever_smoker_no_DBS2 <- ever_smoker %>% filter(DBS2 == 0) #5 

# Ex-smokers and DBS2
ex_smoker_DBS2 <- ex_smoker %>% filter(DBS2>0) #83
ex_smoker_no_DBS2 <- ex_smoker %>% filter(DBS2 == 0) #5 

# Never smokers and DBS2
never_smoker_DBS2 <- never_smoker %>% filter(DBS2>0) #6
never_smoker_no_DBS2 <- never_smoker %>% filter(DBS2 == 0) #2

# Sensitivity and specificity calculations for DBS2
sensitivity_DBS2 <- nrow(ever_smoker_DBS2)/nrow(ever_smoker) #0.9579832
specificity_DBS2 <- nrow(never_smoker_no_DBS2)/nrow(never_smoker) #0.25
PLR_DBS2 <- sensitivity_DBS2/(1-specificity_DBS2) #1.277311
NLR_DBS2 <- (1-sensitivity_DBS2)/specificity_DBS2 #0.1680672
PPV_DBS2 <- nrow(ever_smoker_DBS2)/nrow(DBS2) #0.9344262
NPV_DBS2 <- nrow(never_smoker_no_DBS2)/nrow(no_DBS2) #0.2857143

# ID3 with smoking cohort
# Ever smokers and ID3
ever_smoker_ID3 <- ever_smoker %>% filter(ID3>0) #104
ever_smoker_no_ID3 <- ever_smoker %>% filter(ID3 == 0) #15

# Ex-smokers and ID3
ex_smoker_ID3 <- ex_smoker %>% filter(ID3>0) #74
ex_smoker_no_ID3 <- ex_smoker %>% filter(ID3 == 0) #14 

# Never smokers and ID3
never_smoker_ID3 <- never_smoker %>% filter(ID3>0) #3
never_smoker_no_ID3 <- never_smoker %>% filter(ID3 == 0) #5

# Sensitivity and specificity calculations for ID3
sensitivity_ID3 <- nrow(ever_smoker_ID3)/nrow(ever_smoker) #0.8739496
specificity_ID3 <- nrow(never_smoker_no_ID3)/nrow(never_smoker) #0.625
PLR_ID3 <- sensitivity_ID3/(1-specificity_ID3) #2.330532
NLR_ID3 <- (1-sensitivity_ID3)/specificity_ID3 #0.2016807
PPV_ID3 <- nrow(ever_smoker_ID3)/nrow(ID3) #0.9541284
NPV_ID3 <- nrow(never_smoker_no_ID3)/nrow(no_ID3) #0.25

#Colour palettes for signature plots
my_palette_SBS <- c("black", "grey", "firebrick", "forestgreen", "peachpuff", "orange", "darkturquoise", "lightgoldenrod1","brown",
                           "lightskyblue", "pink", "aquamarine", "orchid1", "lightcoral", "springgreen")

my_palette_DBS <- c("lawngreen", "darkred", "cadetblue1", "pink1", "khaki", "magenta", "yellowgreen")

my_palette_IDS <- c("lightslateblue", "darkseagreen1", "bisque2", "plum2", "slategray2", "tan2", "thistle2","sienna1",  
                     "darkolivegreen2")

# Table 1 - Baseline demographics ------------------------------------------------------
age_summary <- lung_patients_combined %>%
  summarise(
    Median_Age = median(Age.at.diagnosis),
    Age_Range = paste(min(Age.at.diagnosis), "-", max(Age.at.diagnosis))
  )

# Summarize categorical variables
# Gender
sex_summary <- lung_patients_combined %>%
  count(Sex) %>%
  mutate(Proportion = n / sum(n)) %>%
  rename(Category = Sex) %>%
  mutate(Variable = "Sex")

# Stage
# Convert stage to simplified Stage I, II, III, IV
lung_patients_combined <- lung_patients_combined %>%
  mutate(
    Stage = Path.stage..8th.edition., 
    Stage = recode(Stage, 
                   "IA1" = "I", "IA2" = "I", "IA3" = "I", 
                   "IA" = "I", "IB" = "I", "IIA" = "II", 
                   "IIB" = "II", "IIIA" = "III/IV", 
                   "IIIB" = "III/IV", "IVA" = "III/IV")
  )

stage_summary <- lung_patients_combined %>%
  count(Stage) %>%
  mutate(Proportion = n / sum(n)) %>%
  rename(Category = Stage) %>%
  mutate(Variable = "Stage")
  
# Histology
lung_patients_combined <- lung_patients_combined %>%
  mutate(
    Histology = Subtype, 
    Histology = recode(Histology, 
                   "adenosquamous" = "squamous or adenosquamous",
                   "squamous cell" = "squamous or adenosquamous")
  )
histology_summary <- lung_patients_combined %>%
  count(Histology) %>%
  mutate(Proportion = n / sum(n)) %>%
  rename(Category = Histology) %>%
  mutate(Variable = "Histology")

# Performance status (ECOG)
ps_summary <- lung_patients_combined %>%
  count(ECOG.at.diagnosis) %>%
  mutate(Proportion = n / sum(n)) %>%
  rename(Category = ECOG.at.diagnosis) %>%
  mutate(Variable = "PS")
ps_summary$Category <- as.character(ps_summary$Category)

# Smoking status
smoking_status_summary <- lung_patients_combined %>%
  count(Smoking.status) %>%
  mutate(Proportion = n / sum(n)) %>%
  rename(Category = Smoking.status) %>%
  mutate(Variable = "Smoking Status")

# Combine all summaries into one dataframe
final_summary <- bind_rows(
  sex_summary,
  stage_summary,
  histology_summary,
  ps_summary,
  smoking_status_summary
)

# View the result
print(final_summary)

# Figure 1 ---------------------------------------------------------------------
# Figure 1A 
#SBS plot (signature frequency barplot)
# median and range for number of signatures per sample
SBS_sum<-rowSums(lung_patients_combined[,25:39]>0)
summary(SBS_sum)

# plot
SBS_counts <- lung_patients_combined %>%
  summarise(across(starts_with("SBS"), ~ sum(. > 0))) %>%
  pivot_longer(cols = everything(), names_to = "SBS", values_to = "Count")
colnames(SBS_counts) <- c("Signature", "Frequency")

SBS_counts$Signature <- factor(SBS_counts$Signature, levels = SBS_counts$Signature)

SBS_frequency_barplot <- ggplot(SBS_counts, aes(x = Signature, y = Frequency)) + 
  geom_bar(stat = "identity", fill = my_palette_SBS, colour = "black") +
  theme_classic() +
  labs(y = "Patients with signature (n)", x="Mutational signature") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 20), legend.position = "none") +
  theme(axis.title.x = element_text(face = "bold", size = 20)) +
  theme(axis.text.x = element_text(size = 20, angle = 45, vjust = 0.8, hjust = 0.75)) +
  theme(axis.title.y = element_text(face = "bold", size = 20)) +
  theme(axis.text.y = element_text(size = 20)) +
  scale_y_continuous(breaks = seq(0, 140, by = 20), lim = c(0,140))

ggsave2(filename = "Figure_1A.pdf", plot = SBS_frequency_barplot, height = 6, width = 20, units = "in", path = "~/file/path/Figures_final/", device = cairo_pdf)

# Figure 1B
#DBS plot (signature frequency barplot)
# median and range
DBS_sum<-rowSums(lung_patients_combined[,40:46]>0)
summary(DBS_sum)

#plot
DBS_counts <- lung_patients_combined %>%
  summarise(across(starts_with("D"), ~ sum(. > 0))) %>%
  pivot_longer(cols = everything(), names_to = "DBS", values_to = "Count")
colnames(DBS_counts) <- c("Signature", "Frequency")

DBS_counts$Signature <- factor(DBS_counts$Signature, levels = DBS_counts$Signature)

DBS_frequency_barplot <- ggplot(DBS_counts, aes(x = Signature, y = Frequency)) + 
  geom_bar(stat = "identity", fill = my_palette_DBS, colour = "black") +
  theme_classic() +
  labs(y = "Patients with signature (n)", x="Mutational signature") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 20), legend.position = "none") +
  theme(axis.title.x = element_text(face = "bold", size = 20)) +
  theme(axis.text.x = element_text(size = 20, angle = 45, vjust = 0.8, hjust = 0.75)) +
  theme(axis.title.y = element_text(face = "bold", size = 20)) +
  theme(axis.text.y = element_text(size = 20)) +
  scale_y_continuous(breaks = seq(0, 140, by = 20), lim = c(0,140))

ggsave2(filename = "Figure_1B.pdf", plot = DBS_frequency_barplot, height = 6, width = 8, units = "in", path = "~/file/path/Figures_final/", device = cairo_pdf)

#Figure 1C
#ID plot (signature frequency barplot)
# median and range
ID_sum<-rowSums(lung_patients_combined[,47:55]>0)
summary(ID_sum)

# plot
ID_counts <- lung_patients_combined %>%
  summarise(across(starts_with("ID"), ~ sum(. > 0))) %>%
  pivot_longer(cols = everything(), names_to = "ID", values_to = "Count")
colnames(ID_counts) <- c("Signature", "Frequency")

ID_counts$Signature <- factor(ID_counts$Signature, levels = ID_counts$Signature)

ID_frequency_barplot <- ggplot(ID_counts, aes(x = Signature, y = Frequency)) + 
  geom_bar(stat = "identity", fill = my_palette_IDS, colour = "black") +
  theme_classic() +
  labs(y = "Patients with signature (n)", x="Mutational signature") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 20), legend.position = "none") +
  theme(axis.title.x = element_text(face = "bold", size = 20)) +
  theme(axis.text.x = element_text(size = 20, angle = 45, vjust = 0.8, hjust = 0.75)) +
  theme(axis.title.y = element_text(face = "bold", size = 20)) +
  theme(axis.text.y = element_text(size = 20)) +
  scale_y_continuous(breaks = seq(0, 140, by = 20), lim = c(0,140))

ggsave2(filename = "Figure_1C.pdf", plot = ID_frequency_barplot, height = 6, width = 10, units = "in", path = "~/file/path/Figures_final/", device = cairo_pdf)

# Figure 1D - boxplot of SBS mutations per sample where SBS present
#Function to use boxplot.stats to set the box-and-whisker locations  
mybxp = function(x) {
  bxp = log10(boxplot.stats(10^x)[["stats"]])
  names(bxp) = c("ymin","lower", "middle","upper","ymax")
  return(bxp)
}  

# Function to use boxplot.stats for the outliers
myout = function(x) {
  data.frame(y=log10(boxplot.stats(10^x)[["out"]]))
}

COSMIC_SBS_supp <- COSMIC_SBS
COSMIC_SBS_supp_melt <- reshape::melt(COSMIC_SBS_supp)
colnames(COSMIC_SBS_supp_melt) <- c("SBS", "Counts")
COSMIC_SBS_supp_melt <-  COSMIC_SBS_supp_melt %>% filter(Counts>0)

COSMIC_SBS_supp_boxplot <- ggplot(COSMIC_SBS_supp_melt, aes(x = SBS, y = Counts, fill = SBS)) + 
  stat_summary(fun.data=mybxp, geom="boxplot") +
  stat_summary(fun.data=myout, geom="point") +
  theme_classic() +
  labs(y = "Number of mutations", x = "Mutational signature") +
  theme(plot.title = element_blank(), legend.position = "none") +
  theme(axis.title.x = element_text(face = "bold", size = 20)) +
  theme(axis.text.x = element_text(size = 20, angle = 45, vjust = 0.8, hjust = 0.75)) +
  theme(axis.title.y = element_text(face = "bold", size = 20)) +
  theme(axis.text.y = element_text(size = 20)) +
  scale_y_log10(limits = c(1,1000000), breaks = c(1, 10, 100, 1000, 10000, 100000, 1000000, 10000000), labels = function (x) format(x, scientific = F))+
  scale_fill_manual(values= my_palette_SBS) 

ggsave2(filename = "Figure_1D.pdf", plot = COSMIC_SBS_supp_boxplot, height = 6, width = 8, units = "in", path = "~/file/path/Figures_final/", device = cairo_pdf)

# Figure 1E - boxplot of DBS mutations per sample where SBS present
COSMIC_DBS_supp <- COSMIC_DBS
COSMIC_DBS_supp_melt <- reshape::melt(COSMIC_DBS_supp)
colnames(COSMIC_DBS_supp_melt) <- c("DBS", "Counts")
COSMIC_DBS_supp_melt <-  COSMIC_DBS_supp_melt %>% filter(Counts>0)

COSMIC_DBS_supp_boxplot <- ggplot(COSMIC_DBS_supp_melt, aes(x = DBS, y = Counts, fill = DBS)) + 
  stat_summary(fun.data=mybxp, geom="boxplot") +
  stat_summary(fun.data=myout, geom="point") +
  theme_classic() +
  labs(y = "Number of mutations", x = "Mutational signature") +
  theme(plot.title = element_blank(), legend.position = "none") +
  theme(axis.title.x = element_text(face = "bold", size = 20)) +
  theme(axis.text.x = element_text(size = 20, angle = 45, vjust = 0.8, hjust = 0.75)) +
  theme(axis.title.y = element_text(face = "bold", size = 20)) +
  theme(axis.text.y = element_text(size = 20)) +
  scale_y_log10(limits = c(1,1000000), breaks = c(1, 10, 100, 1000, 10000, 100000, 1000000, 10000000), labels = function (x) format(x, scientific = F))+
  scale_fill_manual(values= my_palette_DBS) 

ggsave2(filename = "Figure_1E.pdf", plot = COSMIC_DBS_supp_boxplot, height = 6, width = 8, units = "in", path = "~/file/path/Figures_final/", device = cairo_pdf)

# Figure 1F - boxplot of IDS mutations per sample where SBS present
COSMIC_IDS_supp <- COSMIC_IDS
COSMIC_IDS_supp_melt <- reshape::melt(COSMIC_IDS_supp)
colnames(COSMIC_IDS_supp_melt) <- c("IDS", "Counts")
COSMIC_IDS_supp_melt <-  COSMIC_IDS_supp_melt %>% filter(Counts>0)

COSMIC_IDS_supp_boxplot <- ggplot(COSMIC_IDS_supp_melt, aes(x = IDS, y = Counts, fill = IDS)) + 
  stat_summary(fun.data=mybxp, geom="boxplot") +
  stat_summary(fun.data=myout, geom="point") +
  theme_classic() +
  labs(y = "Number of mutations", x = "Mutational signature") +
  theme(plot.title = element_blank(), legend.position = "none") +
  theme(axis.title.x = element_text(face = "bold", size = 20)) +
  theme(axis.text.x = element_text(size = 20, angle = 45, vjust = 0.8, hjust = 0.75)) +
  theme(axis.title.y = element_text(face = "bold", size = 20)) +
  theme(axis.text.y = element_text(size = 20)) +
  scale_y_log10(limits = c(1,1000000), breaks = c(1, 10, 100, 1000, 10000, 100000, 1000000, 10000000), labels = function (x) format(x, scientific = F))+
  scale_fill_manual(values= my_palette_IDS) 

ggsave2(filename = "Figure_1F.pdf", plot = COSMIC_IDS_supp_boxplot, height = 6, width = 8, units = "in", path = "~/file/path/Figures_final/", device = cairo_pdf)

# Plot Figure 1 A - F as a grid
top_row <- plot_grid(SBS_frequency_barplot, DBS_frequency_barplot, ID_frequency_barplot, nrow = 1)
bottom_row <- plot_grid(COSMIC_SBS_supp_boxplot, COSMIC_DBS_supp_boxplot, COSMIC_IDS_supp_boxplot, nrow = 1)
Figure_1_panel <- plot_grid(top_row, bottom_row, ncol = 1)

ggsave2(filename = "Figure_1_panel.pdf", plot = Figure_1_panel, height = 20, width = 26, units = "in", path = "~/file/path/Figures_final/", device = cairo_pdf)

# Supp Figure 1 ---------------------------------------------------------------------
# Figure S1A TMB whole cohort
# median/range
summary(lung_patients_combined$TMB)

# plot
TMB_melt <- reshape2::melt(lung_patients_combined$TMB)

TMB_cohort_boxplot <- ggplot(TMB_melt, aes(x = "Cohort", y = value, fill = "grey")) + 
  geom_boxplot(fill = "grey", colour = "black", outlier.shape = NA) +
  geom_jitter() +
  theme_classic() +
  labs(x = "", y="TMB (mutations/Mb)") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 20), legend.position = "none") +
  theme(axis.text.x = element_text(size = 20, face = "bold")) +
  theme(axis.title.y = element_text(face = "bold", size = 20)) +
  theme(axis.text.y = element_text(size = 20)) +
  scale_y_continuous(limits = c(0,50), breaks = seq(0, 50, by=10))

ggsave2(filename = "Figure_S1A.pdf", plot = TMB_cohort_boxplot, height = 6, width = 4, units = "in", path = "~/file/path/Figures_final/", device = cairo_pdf)

#Figure S1B TMB squamous vs non-squamous
non_squamous <- lung_patients_combined[which(lung_patients_combined$Subtype == "adenocarcinoma"),]
squamous <- lung_patients_combined[which(lung_patients_combined$Subtype != "adenocarcinoma"),]
lung_patients_combined$Subtype_new <- lung_patients_combined$Subtype
lung_patients_combined$Subtype_new <- gsub("adenocarcinoma", "Non-Squamous", lung_patients_combined$Subtype_new)
lung_patients_combined$Subtype_new <- gsub("adenosquamous", "Squamous", lung_patients_combined$Subtype_new)
lung_patients_combined$Subtype_new <- gsub("squamous cell", "Squamous", lung_patients_combined$Subtype_new)

lung_patients_combined$Subtype_new <- factor(lung_patients_combined$Subtype_new, levels = c("Non-Squamous", "Squamous"))
lung_patients_combined$TMB <- as.numeric(lung_patients_combined$TMB)

summary(non_squamous$TMB)
summary(squamous$TMB)
wilcox.test(non_squamous$TMB, squamous$TMB)

TMB_subtype_stat <- wilcox.test(non_squamous$TMB, squamous$TMB)
my_comparisons <- list(c("Non-Squamous", "Squamous"))

TMB_subtype_boxplot <- ggplot(lung_patients_combined, aes(x=Subtype_new, y=TMB, fill=Subtype_new)) + 
  geom_boxplot(colour = "black", outlier.shape = NA) +
  geom_jitter() +
  theme_classic() +
  labs(x = "", y="TMB (mutations/Mb)") +
  theme(plot.title = element_blank(), legend.position = "none") +
  theme(axis.text.x = element_text(size = 20, face = "bold")) +
  theme(axis.title.y = element_text(face = "bold", size = 20)) +
  theme(axis.text.y = element_text(size = 20)) +
  scale_y_continuous(limits = c(0,55), breaks = seq(0, 55, by=10), expand = c(0,0)) +
  scale_fill_manual(values = c("Non-Squamous" = "purple", "Squamous" = "orange")) +
  geom_signif(comparisons = my_comparisons, textsize=8, fontface = "bold", step_increase = 0.1, map_signif_level = F, 
              vjust = -0.25, annotations = paste("p =", paste(signif(TMB_subtype_stat$p.value, digits = 2))))

ggsave2(filename = "Figure_S1B.pdf", plot = TMB_subtype_boxplot, height = 6, width = 6, units = "in", path = "~/file/path/Figures_final/", device = cairo_pdf)

#Figure S1C
coverage_TMB_stat <- summary(lm(lung_patients_combined$Genome.Wide.Coverage.Mean..X.~lung_patients_combined$TMB))

TMB_coverage_plot <- ggplot(lung_patients_combined, aes(x = Genome.Wide.Coverage.Mean..X., y = TMB))+
  geom_point() +
  theme_classic() +
  labs(x = paste0("Coverage", " (\U00D7)"), y="TMB (mutations/Mb)") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 20), legend.position = "none") +
  theme(axis.title.x = element_text(face = "bold", size = 20)) +
  theme(axis.text.x = element_text(size = 20)) +
  theme(axis.title.y = element_text(face = "bold", size = 20)) +
  theme(axis.text.y = element_text(size = 20)) +
  scale_x_continuous(limits = c(80,160), breaks = seq(80, 160, by=20)) +
  scale_y_continuous(limits = c(0,50), breaks = seq(0, 50, by=10)) +
  geom_smooth(method = lm, se = T, formula = 'y ~ x') + #linear regression 
  annotate("text", x = 140, y = 50, size = 8, label = paste("R\u00b2 adj =", paste(signif(coverage_TMB_stat$adj.r.squared, digits = 2)))) +
  annotate("text", x = 140, y = 47, size = 8, label = paste("p =", paste(signif(coverage_TMB_stat$coefficients[2,4], digits = 2))))

ggsave2(filename = "Figure_S1C.pdf", plot = TMB_coverage_plot, height = 6, width = 6, units = "in", path = "~/file/path/Figures_final/", device = cairo_pdf)

#Figure S1D
purity_TMB_stat <- summary(lm(lung_patients_combined$C_cube.tumour.purity.... ~ lung_patients_combined$TMB))

TMB_purity_plot <- ggplot(lung_patients_combined, aes(x = C_cube.tumour.purity...., y = TMB))+
  geom_point() +
  theme_classic() +
  labs(x = "Tumour purity (%)", y="TMB (mutations/Mb)") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 20), legend.position = "none") +
  theme(axis.title.x = element_text(face = "bold", size = 20)) +
  theme(axis.text.x = element_text(size = 20)) +
  theme(axis.title.y = element_text(face = "bold", size = 20)) +
  theme(axis.text.y = element_text(size = 20)) +
  scale_x_continuous(limits = c(0,100), breaks = seq(0, 100, by = 20)) +
  scale_y_continuous(limits = c(0,50), breaks = seq(0, 50, by=10)) +
  geom_smooth(method = lm, se = T, formula = 'y ~ x') + #linear regression 
  annotate("text", x = 80, y = 50, size = 8, label = paste("R\u00b2 adj =", paste(signif(purity_TMB_stat$adj.r.squared, digits = 2)))) +
  annotate("text", x = 80, y = 47, size = 8, label = paste("p =", paste(signif(purity_TMB_stat$coefficients[2,4], digits = 2))))

ggsave2(filename = "Figure_S1D.pdf", plot = TMB_purity_plot, height = 6, width = 6, units = "in", path = "~/file/path/Figures_final/", device = cairo_pdf)

# Figure S1E - Total number of mutations for each signature across entire cohort
total_SBS <- sum(COSMIC_SBS)
total_DBS <- sum(COSMIC_DBS)
total_IDS <- sum(COSMIC_IDS)

total_mutations <- data.frame(t(data.frame(total_SBS, total_DBS, total_IDS)))
total_mutations$Class <- c("SBS", "DBS", "Indel")
colnames(total_mutations) <- c("Total_mutations", "Mutation_class")
total_mutations$Mutation_class <- factor(total_mutations$Mutation_class, levels = c("SBS", "DBS", "Indel"))

total_mutations_barplot <- ggplot(total_mutations, aes(x = Mutation_class, y = Total_mutations)) + 
  geom_bar(stat = "identity", fill = "grey", colour = "black") +
  theme_classic() +
  labs(y = "Number of mutations", x="Mutation class") +
  theme(plot.title = element_blank(), legend.position = "none") +
  theme(axis.title.x = element_text(face = "bold", size = 20)) +
  theme(axis.text.x = element_text(size = 20)) +
  theme(axis.title.y = element_text(face = "bold", size = 20)) +
  theme(axis.text.y = element_text(size = 20)) +
  scale_y_continuous(breaks = seq(0, 6000000, by = 1000000), lim = c(0,6000000), labels = function (x) format(x, scientific = F))

ggsave2(filename = "Figure_S1E.pdf", plot = total_mutations_barplot, height = 6, width = 4, units = "in", path = "~/file/path/Figures_final/", device = cairo_pdf)

# Figure S1F - total number of mutations for each signature type, for each sample
SBS_per_sample <- data.frame(rowSums(COSMIC_SBS))
SBS_per_sample$Signature <- "SBS"
colnames(SBS_per_sample) <- c("Counts", "Signature")
rownames(SBS_per_sample) <- NULL

DBS_per_sample <- data.frame(rowSums(COSMIC_DBS))
DBS_per_sample$Signature <- "DBS"
colnames(DBS_per_sample) <- c("Counts", "Signature")
rownames(DBS_per_sample) <- NULL

IDS_per_sample <- data.frame(rowSums(COSMIC_IDS))
IDS_per_sample$Signature <- "Indels"
colnames(IDS_per_sample) <- c("Counts", "Signature")
rownames(IDS_per_sample) <- NULL

mutations_per_sample <- rbind(SBS_per_sample, DBS_per_sample, IDS_per_sample)
mutations_per_sample$Signature <- factor(mutations_per_sample$Signature, levels = c("SBS", "DBS", "Indels"))

#Figure S1F
SBS_mutations_per_sample <- ggplot(data = SBS_per_sample, aes(x = Signature, y = Counts)) +
  stat_summary(fun.data=mybxp, geom="boxplot", fill = "grey") +
  stat_summary(fun.data=myout, geom="point", color = "white") +
  geom_jitter() +
  theme_classic() +
  labs(x = "Mutation class", y="Number of mutations") +
  theme(plot.title = element_blank(), legend.position = "none") +
  theme(axis.title.x = element_text(face = "bold", size = 20)) +
  theme(axis.text.x = element_text(size = 20)) +
  theme(axis.title.y = element_text(face = "bold", size = 20)) +
  theme(axis.text.y = element_text(size = 20)) +
  scale_y_log10(limits = c(1,1000000), breaks = c(1, 10, 100, 1000, 10000, 100000, 1000000, 10000000), labels = function (x) format(x, scientific = F))

ggsave2(filename = "Figure_S1F.pdf", plot = SBS_mutations_per_sample, height = 6, width = 4, units = "in", path = "~/file/path/Figures_final/", device = cairo_pdf)

#Figure S1G (log10 transformed)
DBS_mutations_per_sample <- ggplot(data = DBS_per_sample, aes(x = Signature, y = Counts)) +
  stat_summary(fun.data=mybxp, geom="boxplot", fill = "grey") +
  stat_summary(fun.data=myout, geom="point", color = "white") +
  geom_jitter() +
  theme_classic() +
  labs(x = "Mutation class", y="Number of mutations") +
  theme(plot.title = element_blank(), legend.position = "none") +
  theme(axis.title.x = element_text(face = "bold", size = 20)) +
  theme(axis.text.x = element_text(size = 20)) +
  theme(axis.title.y = element_text(face = "bold", size = 20)) +
  theme(axis.text.y = element_text(size = 20)) +
  scale_y_log10(limits = c(1,10000), breaks = c(1, 10, 100, 1000, 10000), labels = function (x) format(x, scientific = F)) 

ggsave2(filename = "Figure_S1G.pdf", plot = DBS_mutations_per_sample, height = 6, width = 4, units = "in", path = "~/file/path/Figures_final/", device = cairo_pdf)

#Figure S1H (log10 transformed)
IDS_mutations_per_sample <- ggplot(data = IDS_per_sample, aes(x = Signature, y = Counts)) +
  stat_summary(fun.data=mybxp, geom="boxplot", fill = "grey") +
  stat_summary(fun.data=myout, geom="point", color = "white") +
  geom_jitter() +
  theme_classic() +
  labs(x = "Mutation Class", y="Number of mutations") +
  theme(plot.title = element_blank(), legend.position = "none") +
  theme(axis.title.x = element_text(face = "bold", size = 20)) +
  theme(axis.text.x = element_text(size = 20)) +
  theme(axis.title.y = element_text(face = "bold", size = 20)) +
  theme(axis.text.y = element_text(size = 20)) +
  scale_y_log10(limits = c(1,100000), breaks = c(1, 10, 100, 1000, 10000, 100000), labels = function (x) format(x, scientific = F)) 

ggsave2(filename = "Figure_S1H.pdf", plot = IDS_mutations_per_sample, height = 6, width = 4, units = "in", path = "~/file/path/Figures_final/", device = cairo_pdf)

# Figure S1 panel
top_row <- plot_grid(TMB_cohort_boxplot, TMB_subtype_boxplot, TMB_coverage_plot, TMB_purity_plot, nrow=1)
bottom_row <- plot_grid(total_mutations_barplot, SBS_mutations_per_sample, DBS_mutations_per_sample, IDS_mutations_per_sample, nrow=1 )
Figure_S1 <- plot_grid(top_row, bottom_row, nrow=2)

ggsave2(filename = "Figure_S1_panel.pdf", plot = Figure_S1, height = 20, width = 26, units = "in", path = "~/file/path/Figures_final/", device = cairo_pdf)

# Figure 2 ----------------------------------------------------------------------

#Barplot of TMB including indels as a percentage, ordered by decreasing SNV count
par(xpd=F, xaxs="i", mfrow = c(1,1), mar=c(4,4,4,4)) 

# TMB barplot
normalised_sigs_TMB <- lung_patients_no_ID[order(lung_patients_no_ID$TMB, decreasing = T),]
normalised_sigs_SBS4 <- lung_patients_no_ID[order(lung_patients_no_ID$SBS4, decreasing = T),]
normalised_sigs_DBS2 <- lung_patients_no_ID[order(lung_patients_no_ID$DBS2, decreasing = T),]
normalised_sigs_ID3 <- lung_patients_no_ID[order(lung_patients_no_ID$ID3, decreasing = T),]

TMB_plot <- ggplot(normalised_sigs_TMB, aes(x = reorder(Gel_ID, -TMB), y = TMB)) +
  geom_bar(stat = "identity", show.legend = F, colour = "black") +
  theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  labs(y = "TMB (Mutations/Mb)") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 20)) +
  theme(axis.text.x = element_blank()) +
  theme(axis.title.x = element_blank()) +
  theme(axis.title.y = element_text(face = "bold", size = 14, vjust = 0.5)) +
  theme(axis.text.y = element_text(size = 10)) +
  theme(axis.ticks.x = element_blank()) +
  scale_y_continuous(limits = c(0,50), breaks = seq(0, 50, by = 10), expand = c(0,0))

# Normalised SBS plot
SBS_normalised_signatures_df <- SBS[match(normalised_sigs_TMB$Gel_ID, SBS$Gel_ID),] # here need to change 
SBS_normalised_signatures_df[is.na(SBS_normalised_signatures_df)] <- 0
rownames(SBS_normalised_signatures_df) <- NULL
SBS_normalised_signatures_df$Gel_ID <- factor(SBS_normalised_signatures_df$Gel_ID, levels = SBS_normalised_signatures_df$Gel_ID)
SBS_normalised_signatures_df_melt <- reshape2::melt(SBS_normalised_signatures_df)
SBS_normalised_signatures_df_melt$value <- 100*SBS_normalised_signatures_df_melt$value

SBS_TMB_plot <- ggplot(SBS_normalised_signatures_df_melt, aes(x = Gel_ID, y = value, fill = variable, group = value)) +
  geom_bar(stat = "identity", show.legend = T,colour = "black") +
  theme_classic() + 
  scale_fill_manual(values = my_palette_SBS) +
  labs(y = "SBS (%)", x = "Patient") +
  theme(plot.title = element_blank()) +
  theme(axis.title.x = element_blank()) +
  theme(axis.text.x = element_blank()) +
  theme(axis.title.y = element_text(face = "bold", size = 14, vjust = 0.5)) +
  theme(axis.text.y = element_text(size = 10)) +
  scale_y_continuous(breaks = seq(0, 100, by = 20), expand = c(0,0)) +
  theme(legend.position = "none")

# Normalised DBS plot
DBS_normalised_signatures_df <- DBS[match(normalised_sigs_TMB$Gel_ID, DBS$Gel_ID),]
DBS_normalised_signatures_df[is.na(DBS_normalised_signatures_df)] <- 0
rownames(DBS_normalised_signatures_df) <- NULL
DBS_normalised_signatures_df$Gel_ID <- factor(DBS_normalised_signatures_df$Gel_ID, levels = DBS_normalised_signatures_df$Gel_ID)
DBS_normalised_signatures_df_melt <- reshape2::melt(DBS_normalised_signatures_df)
DBS_normalised_signatures_df_melt$value <- 100*DBS_normalised_signatures_df_melt$value

DBS_TMB_plot <- ggplot(DBS_normalised_signatures_df_melt, aes(x = Gel_ID, y = value, fill = variable, group = value)) +
  geom_bar(stat = "identity", show.legend = T, colour = "black") +
  theme_classic() + 
  scale_fill_manual(values = my_palette_DBS) +
  labs(y = "DBS (%)", x = "Patient") +
  theme(plot.title = element_blank()) +
  theme(axis.title.x = element_blank()) +
  theme(axis.text.x = element_blank()) +
  theme(axis.title.y = element_text(face = "bold", size = 14, vjust = 0.5)) +
  theme(axis.text.y = element_text(size = 10)) +
  scale_y_continuous(breaks = seq(0, 100, by = 20), expand = c(0,0)) +
  theme(legend.position = "none")

# Normalised IDS plot
IDS_normalised_signatures_df <- IDS[match(normalised_sigs_TMB$Gel_ID, IDS$Gel_ID),]
IDS_normalised_signatures_df[is.na(IDS_normalised_signatures_df)] <- 0
rownames(IDS_normalised_signatures_df) <- NULL
IDS_normalised_signatures_df$Gel_ID <- factor(IDS_normalised_signatures_df$Gel_ID, levels = IDS_normalised_signatures_df$Gel_ID)
IDS_normalised_signatures_df_melt <- reshape2::melt(IDS_normalised_signatures_df)
IDS_normalised_signatures_df_melt$value <- 100*IDS_normalised_signatures_df_melt$value

IDS_TMB_plot <- ggplot(IDS_normalised_signatures_df_melt, aes(x = Gel_ID, y = value, fill = variable, group = value)) +
  geom_bar(stat = "identity", show.legend = T, colour = "black") +
  theme_classic() + 
  scale_fill_manual(values = my_palette_IDS) +
  labs(y = "Indel (%)", x = "Patient") +
  theme(plot.title = element_blank()) +
  theme(axis.title.x = element_blank()) +
  theme(axis.text.x = element_blank()) +
  theme(axis.title.y = element_text(face = "bold", size = 14, vjust = 0.5)) +
  theme(axis.text.y = element_text(size = 10)) +
  scale_y_continuous(breaks = seq(0, 100, by = 20), expand = c(0,0)) +
  theme(legend.position = "none")

#Smoking status (one line thick)
smoking_heatmap <- normalised_sigs_TMB[,c("Gel_ID", "Smoking.status", "Smoking.status"),]
colnames(smoking_heatmap) <- c("Gel_ID", "Smoking_status", "Colour")
smoking_heatmap$Smoking_status[smoking_heatmap$Smoking_status=="Current"] <- 1
smoking_heatmap$Smoking_status[smoking_heatmap$Smoking_status=="Ex"] <- 2
smoking_heatmap$Smoking_status[smoking_heatmap$Smoking_status=="Never"] <- 3

smoking_heatmap$Colour[smoking_heatmap$Colour=="Current"] <- "black"
smoking_heatmap$Colour[smoking_heatmap$Colour=="Ex"] <- "grey"
smoking_heatmap$Colour[smoking_heatmap$Colour=="Never"] <- "pink"
smoking_heatmap$Colour[is.na(smoking_heatmap$Colour)] <- "gray88"
        
smoking_heatmap$Smoking_status <- 1
smoking_heatmap$Gel_ID <- factor(smoking_heatmap$Gel_ID, levels = smoking_heatmap$Gel_ID)
      
smoking_ggplot <- ggplot(data = smoking_heatmap, aes(x = Gel_ID, y = Smoking_status, fill = Colour)) +
  geom_bar(stat = "identity", show.legend = T, colour = "black") +
  scale_y_continuous(expand = c(0,0)) +
  theme_void()+
  scale_fill_manual(label = c("Current", "Ex", "Never", "Unknown"), values = c("black" = "black", "grey" = "grey", "pink" = "pink", "gray88" = "gray88")) +
  labs(y = "Smoking status") +
  theme(axis.title.y = element_text(face = "bold", size = 14, angle = 0, vjust = 0.5)) +
  theme(legend.position = "none") 
      
#Gender plot
sex_heatmap <- normalised_sigs_TMB[,c("Gel_ID", "Sex", "Sex"),]
colnames(sex_heatmap) <- c("Gel_ID", "Sex", "Colour")
sex_heatmap$Sex[sex_heatmap$Sex=="M"] <- 1
sex_heatmap$Sex[sex_heatmap$Sex=="F"] <- 2
      
sex_heatmap$Colour[sex_heatmap$Colour=="M"] <- "firebrick"
sex_heatmap$Colour[sex_heatmap$Colour=="F"] <- "steelblue"
          
sex_heatmap$Sex <- 1
sex_heatmap$Gel_ID <- factor(sex_heatmap$Gel_ID, levels = sex_heatmap$Gel_ID)
        
sex_ggplot <- ggplot(data = sex_heatmap, aes(x = Gel_ID, y = Sex, fill = Colour)) +
   geom_bar(stat = "identity", show.legend = T, colour = "black") +
   scale_y_continuous(expand = c(0,0)) +
   theme_void()+
   scale_fill_manual(labels = c("Male", "Female"), values = c("firebrick" = "firebrick", "steelblue" = "steelblue")) +
   labs(y = "Sex") +
   theme(axis.title.y = element_text(face = "bold", size = 14, angle = 0, vjust = 0.5)) +
   theme(legend.position = "none") 
        
#Subtype status (one line thick)
subtype_heatmap <- normalised_sigs_TMB[,c("Gel_ID", "Subtype", "Subtype"),]
colnames(subtype_heatmap) <- c("Gel_ID", "Subtype", "Colour")
subtype_heatmap$Subtype[subtype_heatmap$Subtype=="adenocarcinoma"] <- 1
subtype_heatmap$Subtype[subtype_heatmap$Subtype=="squamous cell"] <- 2
subtype_heatmap$Subtype[subtype_heatmap$Subtype=="adenosquamous"] <- 2
        
subtype_heatmap$Colour[subtype_heatmap$Colour=="adenocarcinoma"] <- "purple"
subtype_heatmap$Colour[subtype_heatmap$Colour=="squamous cell"] <- "orange"
subtype_heatmap$Colour[subtype_heatmap$Colour=="adenosquamous"] <- "orange"
              
subtype_heatmap$Subtype <- 1
subtype_heatmap$Gel_ID <- factor(subtype_heatmap$Gel_ID, levels = subtype_heatmap$Gel_ID)
            
subtype_ggplot <- ggplot(data = subtype_heatmap, aes(x = Gel_ID, y = Subtype, fill = Colour)) +
    geom_bar(stat = "identity", show.legend = T, colour = "black") +
    scale_y_continuous(expand = c(0,0)) +
    theme_void()+
    scale_fill_manual(labels = c("Non-squamous", "Squamous cell"), values = c("purple" = "purple", "orange" = "orange")) +
    labs(y = "Subtype") +
    theme(axis.title.y = element_text(face = "bold", size = 14, angle = 0, vjust = 0.5))  +
    theme(legend.position = "none") 
            
#Pack_year status (one line thick)
py_heatmap <- normalised_sigs_TMB[,c("Gel_ID", "Smoking..pack.years.", "Smoking..pack.years."),]
colnames(py_heatmap) <- c("Gel_ID", "Pack_years", "Colour")
            
normalised_sigs_TMB$Smoking..pack.years.
py_heatmap$Colour <- py_heatmap$Pack_years
            
py_heatmap$Colour[is.na(py_heatmap$Pack_years)] <- "cornsilk"
py_heatmap$Colour[py_heatmap$Pack_years ==0] <- "blanchedalmond"
py_heatmap$Colour[py_heatmap$Pack_years >0 & py_heatmap$Pack_years<=20] <- "darkseagreen1"
py_heatmap$Colour[py_heatmap$Pack_years > 20 & py_heatmap$Pack_years<=40] <- "darkseagreen2"
py_heatmap$Colour[py_heatmap$Pack_years > 40 & py_heatmap$Pack_years<= 60] <- "darkseagreen3"
py_heatmap$Colour[py_heatmap$Pack_years > 60] <- "darkseagreen"
                        
py_heatmap$Pack_years <- 1
py_heatmap$Gel_ID <- factor(py_heatmap$Gel_ID, levels = py_heatmap$Gel_ID)
                      
py_ggplot <- ggplot(data = py_heatmap, aes(x = Gel_ID, y = Pack_years, fill = Colour)) +
     geom_bar(stat = "identity", show.legend = T, colour = "black") +
     scale_y_continuous(expand = c(0,0)) +
     theme_void()+
     scale_fill_manual(labels = c("Missing", "0py", paste0("\u226420", "py"), "20-40py", "40-60py", paste0("\u02C3 60", "py")), 
                                          values = c("cornsilk" = "cornsilk", "blanchedalmond" = "blanchedalmond", "darkseagreen1" = "darkseagreen1", "darkseagreen2" = "darkseagreen2", "darkseagreen3" = "darkseagreen3", "darkseagreen" = "darkseagreen")) +
     labs(y = "Pack years") +
     theme(axis.title.y = element_text(face = "bold", size = 14, angle = 0, vjust = 0.5))  +
     theme(legend.position = "none") 
                      
x_axis_only_plot <- ggplot(IDS_normalised_signatures_df_melt, aes(x = Gel_ID, y = value, fill = variable, group = value)) +
     geom_blank() +
     theme_bw() +
     theme(axis.line.y=element_blank(),
           axis.text.y=element_blank(),
           axis.ticks.y=element_blank(),
           axis.title.y=element_blank(),
           panel.grid.minor.y=element_blank(),
           panel.grid.minor.x=element_blank(),
           panel.grid.major.x=element_blank(), 
           panel.grid.major.y=element_blank(), 
           panel.border = element_blank(),
           axis.text.x = element_text(size = 10, angle = 90, vjust = 0.5, hjust = 0.95),
           axis.title.x = element_blank(),
           axis.line = element_line("black"),
           axis.ticks.x = element_line("black"),
           plot.margin = unit(c(0, 0, 2,0), "cm"))
                      
SBS_legend <- get_legend(SBS_TMB_plot + 
                           theme(legend.position = "bottom", legend.direction = "horizontal", 
                                 legend.box = "horizontal", legend.title = element_blank(),
                                 legend.text = element_text(size=14),
                                 plot.margin = unit(c(0,0,0,0), "cm")))
as_ggplot(SBS_legend)
                      
DBS_legend <- get_legend(DBS_TMB_plot + theme(legend.position = "bottom", legend.direction = "horizontal", 
                                              legend.box = "horizontal", legend.title = element_blank(),
                                              legend.text = element_text(size=14),
                                              plot.margin = unit(c(0,0,0,0), "cm")))
as_ggplot(DBS_legend)
                      
IDS_legend <- get_legend(IDS_TMB_plot + theme(legend.position = "bottom", legend.direction = "horizontal", 
                                              legend.box = "horizontal", legend.title = element_blank(),
                                              legend.text = element_text(size=14),
                                              plot.margin = unit(c(0,0,0,0), "cm")))
as_ggplot(IDS_legend)
                      
sex_legend <- get_legend(sex_ggplot + theme(legend.position = "bottom", legend.direction = "vertical", legend.box = "vertical",
                                            legend.title = element_blank(), legend.text = element_text(size=14), plot.margin = unit(c(0,0,0,0), "cm")))
as_ggplot(sex_legend)

subtype_legend <- get_legend(subtype_ggplot + theme(legend.position = "bottom", legend.direction = "vertical", 
                                                    legend.box = "vertical", legend.title = element_blank(),
                                                    legend.text = element_text(size=14),
                                                    plot.margin = unit(c(0,0,0,0), "cm")))
as_ggplot(subtype_legend)

smoking_legend <- get_legend(smoking_ggplot + theme(legend.position = "bottom", legend.direction = "vertical", 
                                                    legend.box = "vertical", legend.title = element_blank(), 
                                                    legend.text = element_text(size=14),
                                                    plot.margin = unit(c(0,0,0,0), "cm")))
as_ggplot(smoking_legend)

py_legend <- get_legend(py_ggplot + theme(legend.position = "bottom", legend.direction = "vertical", 
                                          legend.box = "vertical", legend.title = element_blank(), 
                                          legend.text = element_text(size=14),
                                          plot.margin = unit(c(0,0,0,0), "cm")))
as_ggplot(py_legend)

legends_1 <- plot_grid(SBS_legend, DBS_legend, IDS_legend, align = "h", nrow = 3)
legends_2 <- plot_grid(sex_legend, subtype_legend, smoking_legend, py_legend, align = "h", nrow = 1)

all_plots <- plot_grid(TMB_plot,SBS_TMB_plot, DBS_TMB_plot, IDS_TMB_plot, sex_ggplot, subtype_ggplot, smoking_ggplot, py_ggplot,
                       x_axis_only_plot, align = "v", ncol = 1, rel_heights = c(1,1,1,1,0.1,0.1,0.1,0.1,0.6))


ggsave2(filename = "Figure_2_SBS_legend2.pdf", plot = legends_1, height = 12, width = 15, units = "in", path = "~/file/path/Figures_final/", device = cairo_pdf)
ggsave2(filename = "Figure_2_Clinical_legend2.pdf", plot = legends_2, height = 12, width = 15, units = "in", path = "~/file/path/Figures_final/", device = cairo_pdf)
ggsave2(filename = "Figure_2_Main2.pdf", plot = all_plots, height = 14, width = 15, units = "in", path = "~/file/path/Figures_final/", device = cairo_pdf)


# Supp Figure 2 ----------------------------------------------------------------------
# Supp Fig 2A 
SBS4_TMB <- cbind("SBS4-positive", SBS4$TMB)
no_SBS4_TMB <- cbind("SBS4-negative", no_SBS4$TMB)

#DBS2 and TMB boxplots
DBS2_TMB <- cbind("DBS2-positive", DBS2$TMB)
no_DBS2_TMB <- cbind("DBS2-negative", no_DBS2$TMB)

#ID3 and TMB boxplots
ID3_TMB <- cbind("ID3-positive", ID3$TMB)
no_ID3_TMB <- cbind("ID3-negative", no_ID3$TMB)

#Combined smoking signatures and and TMB boxplots
smoking_sigs_TMB_boxplots <- data.frame(rbind(SBS4_TMB, no_SBS4_TMB, DBS2_TMB, no_DBS2_TMB, ID3_TMB, no_ID3_TMB))
colnames(smoking_sigs_TMB_boxplots) <- c("Signatures", "TMB")
smoking_sigs_TMB_boxplots$TMB <- as.numeric(smoking_sigs_TMB_boxplots$TMB)
smoking_sigs_TMB_boxplots$Signatures <- factor(smoking_sigs_TMB_boxplots$Signatures, 
                                               levels = c("SBS4-positive", "SBS4-negative", "DBS2-positive", 
                                                          "DBS2-negative", "ID3-positive", "ID3-negative"))

my_comparisons <- split(t(combn(levels(smoking_sigs_TMB_boxplots$Signatures), 2)), seq(nrow(t(combn(levels(smoking_sigs_TMB_boxplots$Signatures), 2))))) #all pairwise combinations
my_comparisons <- my_comparisons[c(1,10,15)] #pairwise combinations manual order

TMB_smoking_signature_stat <- pairwise.wilcox.test(smoking_sigs_TMB_boxplots$TMB, smoking_sigs_TMB_boxplots$Signatures, p.adjust.method = "BH")

smoking_sigs_TMB_boxplots_ggplot <- ggplot(smoking_sigs_TMB_boxplots, aes(x=Signatures, y=TMB, fill=Signatures)) + 
  geom_boxplot(colour = "black",  outlier.shape = NA) +
  geom_jitter() +
  theme_classic() +
  labs(x = "Mutational signature cohort", y="TMB (mutations/Mb)") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 20), legend.position = "none") +
  theme(axis.title.x = element_text(face = "bold", size = 20)) +
  theme(axis.text.x = element_text(size = 20)) +
  theme(axis.title.y = element_text(face = "bold", size = 20)) +
  theme(axis.text.y = element_text(size = 20)) +
  scale_y_continuous(limits = c(0,55), breaks = seq(0, 55, by = 10), expand = c(0,0)) +
  scale_fill_manual(values = c("firebrick", "grey", "darkred", "grey", "bisque2", "grey")) +
  geom_signif(comparisons = my_comparisons, textsize=8, fontface = "bold", y_position = 47, vjust = -0.5, 
              annotations = c(paste("p =", paste(signif(TMB_smoking_signature_stat$p.value[1,1], digits = 2))),
                              paste("p =", paste(signif(TMB_smoking_signature_stat$p.value[3,3], digits = 2))),
                              paste("p =", paste(signif(TMB_smoking_signature_stat$p.value[5,5], digits = 2))))) #annotations means manually enter the p-values

ggsave2(filename = "Figure_S2A.pdf", plot = smoking_sigs_TMB_boxplots_ggplot, height = 6, width = 15, units = "in", path = "~/file/path/Figures_final/", device = cairo_pdf)


# Figure 3 ------------------------------#####
# Figure 3A 
# Venn diagrams
venn_df <- lung_patients %>%
  select(c("SBS4", "DBS2", "ID3", "Smoking.status"))
venn_df[venn_df==0] <- NA

# Figure 3A - entire cohort

area1 <- nrow(unique(venn_df[which(!is.na(venn_df$SBS4) & is.na(venn_df$DBS2) 
                                   & is.na(venn_df$ID3)),]))
area2 <- nrow(unique(venn_df[which(is.na(venn_df$SBS4) & !is.na(venn_df$DBS2) 
                                   & is.na(venn_df$ID3)),]))
area3 <- nrow(unique(venn_df[which(is.na(venn_df$SBS4) & is.na(venn_df$DBS2) 
                                   & !is.na(venn_df$ID3)),]))

n12 <-nrow(unique(venn_df[which(!is.na(venn_df$SBS4) & !is.na(venn_df$DBS2) 
                                & is.na(venn_df$ID3)),]))
n13 <- nrow(unique(venn_df[which(!is.na(venn_df$SBS4) & is.na(venn_df$DBS2) 
                                 & !is.na(venn_df$ID3)),]))
n23 <- nrow(unique(venn_df[which(is.na(venn_df$SBS4) & !is.na(venn_df$DBS2) 
                                 & !is.na(venn_df$ID3)),]))

n123 <- nrow(unique(venn_df[which(!is.na(venn_df$SBS4) & !is.na(venn_df$DBS2) 
                                  & !is.na(venn_df$ID3)),]))

fit1 <- c("area1" = area1, "area2" = area2, "area3" = area3, "area1&area2" = n12, "area1&area3" = n13, 
          "area2&area3" = n23, "area1&area2&area3" = n123)

venn <- venn(fit1)
eulerr_options(labels = list(cex = 2))
euler_1 <- plot(venn, labels = c("SBS4", "DBS2", "ID3"), quantities = list(cex = 2), 
                fills = c("grey", "peachpuff", "lightblue"), fontface = "bold")

ggsave2(filename = "Figure_3A_all.pdf", plot = euler_1, height = 6, width = 6, units = "in", path = "~/file/path/Figures_final/", device = cairo_pdf)

# Figure 3A ever smoker 
area1 <- nrow(unique(venn_df[which(!is.na(venn_df$SBS4) & is.na(venn_df$DBS2) 
                                   & is.na(venn_df$ID3) & venn_df$Smoking.status %in% c("Current", "Ex")),]))
area2 <- nrow(unique(venn_df[which(is.na(venn_df$SBS4) & !is.na(venn_df$DBS2) 
                                   & is.na(venn_df$ID3) & venn_df$Smoking.status %in% c("Current", "Ex")),]))
area3 <- nrow(unique(venn_df[which(is.na(venn_df$SBS4) & is.na(venn_df$DBS2) 
                                   & !is.na(venn_df$ID3) & venn_df$Smoking.status %in% c("Current", "Ex")),]))

n12 <-nrow(unique(venn_df[which(!is.na(venn_df$SBS4) & !is.na(venn_df$DBS2) 
                                & is.na(venn_df$ID3) & venn_df$Smoking.status %in% c("Current", "Ex")),]))
n13 <- nrow(unique(venn_df[which(!is.na(venn_df$SBS4) & is.na(venn_df$DBS2) 
                                 & !is.na(venn_df$ID3) & venn_df$Smoking.status %in% c("Current", "Ex")),]))
n23 <- nrow(unique(venn_df[which(is.na(venn_df$SBS4) & !is.na(venn_df$DBS2) 
                                 & !is.na(venn_df$ID3) & venn_df$Smoking.status %in% c("Current", "Ex")),]))

n123 <- nrow(unique(venn_df[which(!is.na(venn_df$SBS4) & !is.na(venn_df$DBS2) 
                                  & !is.na(venn_df$ID3) & venn_df$Smoking.status %in% c("Current", "Ex")),]))

fit2 <- c("area1" = area1, "area2" = area2, "area3" = area3, "area1&area2" = n12, "area1&area3" = n13, 
          "area2&area3" = n23, "area1&area2&area3" = n123)

venn2 <- venn(fit2)
euler_2 <- plot(venn2, labels = c("SBS4", "DBS2", "ID3"), quantities = list(cex = 2), fills = c("grey", "peachpuff", "lightblue"))

ggsave2(filename = "Figure_3A_ever.pdf", plot = euler_2, height = 6, width = 6, units = "in", path = "~/file/path/Figures_final/", device = cairo_pdf)

# Figure 3A Never
area1 <- nrow(unique(venn_df[which(!is.na(venn_df$SBS4) & is.na(venn_df$DBS2) 
                                   & is.na(venn_df$ID3) & venn_df$Smoking.status == "Never"),]))

area2 <- nrow(unique(venn_df[which(is.na(venn_df$SBS4) & !is.na(venn_df$DBS2) 
                                   & is.na(venn_df$ID3) & venn_df$Smoking.status == "Never"),]))
area3 <- nrow(unique(venn_df[which(is.na(venn_df$SBS4) & is.na(venn_df$DBS2) 
                                   & !is.na(venn_df$ID3) & venn_df$Smoking.status == "Never"),]))

n12 <-nrow(unique(venn_df[which(!is.na(venn_df$SBS4) & !is.na(venn_df$DBS2) 
                                & is.na(venn_df$ID3) & venn_df$Smoking.status == "Never"),]))
n13 <- nrow(unique(venn_df[which(!is.na(venn_df$SBS4) & is.na(venn_df$DBS2) 
                                 & !is.na(venn_df$ID3) & venn_df$Smoking.status == "Never"),]))
n23 <- nrow(unique(venn_df[which(is.na(venn_df$SBS4) & !is.na(venn_df$DBS2) 
                                 & !is.na(venn_df$ID3) & venn_df$Smoking.status == "Never"),]))

n123 <- nrow(unique(venn_df[which(!is.na(venn_df$SBS4) & !is.na(venn_df$DBS2) 
                                  & !is.na(venn_df$ID3) & venn_df$Smoking.status == "Never"),]))

fit3 <- c("area1" = area1, "area2" = area2, "area3" = area3, "area1&area2" = n12, "area1&area3" = n13, 
          "area2&area3" = n23, "area1&area2&area3" = n123)

venn3 <- venn(fit3)
euler_3 <- plot(venn3, labels = c("SBS4", "DBS2", "ID3"), quantities = list(cex = 2), fills = c("grey", "peachpuff", "lightblue"))

ggsave2(filename = "Figure_3A_never.pdf", plot = euler_3, height = 6, width = 6, units = "in", path = "~/file/path/Figures_final/", device = cairo_pdf)

# Venn combined
euler_combined <- plot_grid(euler_1, euler_2,euler_3, align = "h", nrow = 1)
ggsave2(filename = "Figure_3A.pdf", plot = euler_combined, height = 6, width = 18, units = "in", path = "~/file/path/Figures_final/", device = cairo_pdf)

# Figure 3B
# Create table with summary sensitivity and specificity calculations 
# This is drawn from calculations at the top of the code (lines 128-195)
results_df <- data.frame(
  Sensitivity = round(c(sensitivity_SBS4, sensitivity_DBS2, sensitivity_ID3),4) *100,
  Specificity = round(c(specificity_SBS4, specificity_DBS2, specificity_ID3),4) *100,
  PPV = round(c(PPV_SBS4, PPV_DBS2, PPV_ID3),4) *100,
  NPV = round(c(NPV_SBS4, NPV_DBS2, NPV_ID3),4) *100,
  PLR = round(c(PLR_SBS4, PLR_DBS2, PLR_ID3),2),
  NLR = round(c(NLR_SBS4, NLR_DBS2, NLR_ID3),2)
)
rownames(results_df)<- c("SBS4", "DBS2", "ID3")
colnames(results_df)[colnames(results_df) == "Specificity"] <- "Specificity (%)"
colnames(results_df)[colnames(results_df) == "Sensitivity"] <- "Sensitivity (%)"
colnames(results_df)[colnames(results_df) == "PPV"] <- "PPV (%)"
colnames(results_df)[colnames(results_df) == "NPV"] <- "NPV (%)"
table_grob <- tableGrob(results_df)
pdf("Figure_3B.pdf", width = 6, height = 2)  # Set PDF dimensions (adjust as needed)
grid.draw(table_grob)  # Draw the table on the PDF device
dev.off()  # Close the PDF device

# Figure 3C
#TMB and smoking status
summary(current_smoker$TMB)
summary(ex_smoker$TMB)
summary(never_smoker$TMB)
summary(ever_smoker$TMB)

current_smoker_df <- cbind("Current", current_smoker$TMB)
ex_smoker_df <- cbind("Ex", ex_smoker$TMB)
never_smoker_df <- cbind("Never", never_smoker$TMB)

TMB_smoking_status <- data.frame(rbind(current_smoker_df, ex_smoker_df, never_smoker_df))
colnames(TMB_smoking_status) <- c("Smoking_status", "TMB")
TMB_smoking_status$Smoking_status <- as.factor(TMB_smoking_status$Smoking_status)
TMB_smoking_status$TMB <- as.numeric(TMB_smoking_status$TMB)
TMB_smoking_status$Smoking_status <- factor(TMB_smoking_status$Smoking_status, levels = c("Current", "Ex", "Never"))

my_comparisons <- split(t(combn(levels(TMB_smoking_status$Smoking_status), 2)), c(1,3,2)) #all pairwise combinations manual order

TMB_smoking_status_stat <- pairwise.wilcox.test(TMB_smoking_status$TMB, TMB_smoking_status$Smoking_status, p.adjust.method = "BH")

TMB_smoking_status_boxplot <- ggplot(TMB_smoking_status, aes(x=Smoking_status, y=TMB)) + 
  geom_boxplot(outlier.shape = NA, fill = "lightgrey") +
  geom_jitter(aes(colour = Smoking_status), size = 3) +
  theme_classic() +
  labs(x = "Smoking status", y="TMB (mutations/Mb)") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 20), 
        legend.title = element_text(face = "bold", size = 15),legend.text = element_text(size = 15),
        axis.title.x = element_text(face = "bold", size = 20), axis.text.x = element_text(size = 15), 
        axis.title.y = element_text(face = "bold", size = 20), axis.text.y = element_text(size = 15)) +
  scale_y_continuous(limits = c(0,60), breaks = seq(0, 60, by=5)) +
  scale_colour_manual(values = c("Current" = "firebrick", "Ex" = "steelblue", 
                                 "Never" = "forestgreen"), name = "Smoking status") +
  geom_signif(comparisons = my_comparisons, textsize=6, fontface = "bold", step_increase = 0.1, map_signif_level = F,
              annotations = c(paste("p =", paste(signif(TMB_smoking_status_stat$p.value[1], digits = 2))),
                              paste("p =", paste(signif(TMB_smoking_status_stat$p.value[4], digits = 2))),
                              paste("p =", paste(signif(TMB_smoking_status_stat$p.value[2], digits = 2))))) 

ggsave2(filename = "Figure_3C.pdf", plot = TMB_smoking_status_boxplot, height = 8, width = 8, units = "in", path =  "~/file/path/Figures_final/", device = cairo_pdf)

#Figure 3D 
#SBS4 and smoking status  

# SBS4 raw counts 
df <- COSMIC_SBS %>%
  tibble::rownames_to_column(var="Gel_ID") 
df$Gel_ID <- as.numeric(df$Gel_ID)

current_smoker_SBS4_activity <- df %>%
  merge(current_smoker, by = "Gel_ID") %>%
  select(Smoking.status, SBS4.x) %>%
  rename(SBS4 = SBS4.x)

ex_smoker_SBS4_activity <- df %>%
  merge(ex_smoker, by = "Gel_ID") %>%
  select(Smoking.status, SBS4.x) %>%
  rename(SBS4 = SBS4.x)

never_smoker_SBS4_activity <- df %>%
  merge(never_smoker, by = "Gel_ID") %>%
  select(Smoking.status, SBS4.x) %>%
  rename(SBS4 = SBS4.x)

SBS4_activity_smoking_status <- data.frame(rbind(current_smoker_SBS4_activity, ex_smoker_SBS4_activity, never_smoker_SBS4_activity))
colnames(SBS4_activity_smoking_status) <- c("Smoking_status", "SBS4")
SBS4_activity_smoking_status$Smoking_status <- as.factor(SBS4_activity_smoking_status$Smoking_status)
SBS4_activity_smoking_status$SBS4 <- as.numeric(SBS4_activity_smoking_status$SBS4)
SBS4_activity_smoking_status$Smoking_status <- factor(SBS4_activity_smoking_status$Smoking_status, levels = c("Current", "Ex", "Never"))

my_comparisons <- split(t(combn(levels(SBS4_activity_smoking_status$Smoking_status), 2)), c(1,3,2)) #all pairwise combinations manual order

SBS4_smoking_status_stat <- pairwise.wilcox.test(SBS4_activity_smoking_status$SBS4, SBS4_activity_smoking_status$Smoking_status, p.adjust.method = "BH")

SBS4_smoking_status_boxplot <- ggplot(SBS4_activity_smoking_status, aes(x=Smoking_status, y=SBS4)) + 
  geom_boxplot(outlier.shape = NA, fill = "lightgrey") +
  geom_jitter(aes(colour = Smoking_status), size = 3) +
  theme_classic() +
  labs(x = "Smoking status", y="SBS4 mutations") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 20), 
        legend.title = element_text(face = "bold", size = 15),legend.text = element_text(size = 15),
        axis.title.x = element_text(face = "bold", size = 20), axis.text.x = element_text(size = 15), 
        axis.title.y = element_text(face = "bold", size = 20), axis.text.y = element_text(size = 15)) +
  scale_y_continuous(limits = c(0,200000), breaks = seq(0, 200000, by=20000)) +
  scale_colour_manual(values = c("Current" = "firebrick", "Ex" = "steelblue", 
                                 "Never" = "forestgreen"), name = "Smoking status") +
  geom_signif(comparisons = my_comparisons, textsize=6, fontface = "bold", step_increase = 0.1, map_signif_level = F,
              annotations = c(paste("p =", paste(signif(SBS4_smoking_status_stat$p.value[1], digits = 2))),
                              paste("p =", paste(signif(SBS4_smoking_status_stat$p.value[4], digits = 2))),
                              paste("p =", paste(signif(SBS4_smoking_status_stat$p.value[2], digits = 2))))) 

ggsave2(filename = "Figure_3D.pdf", plot = SBS4_smoking_status_boxplot, height = 8, width = 8, units = "in", path =  "~/file/path/Figures_final/", device = cairo_pdf)

# Fig 3E
# DBS2
df2 <- COSMIC_DBS %>%
  tibble::rownames_to_column(var="Gel_ID") 
df2$Gel_ID <- as.numeric(df2$Gel_ID)

current_smoker_DBS2_activity <- df2 %>%
  merge(current_smoker, by = "Gel_ID") %>%
  select(Smoking.status, DBS2.x) %>%
  rename(DBS2 = DBS2.x)

ex_smoker_DBS2_activity <- df2 %>%
  merge(ex_smoker, by = "Gel_ID") %>%
  select(Smoking.status, DBS2.x) %>%
  rename(DBS2 = DBS2.x)

never_smoker_DBS2_activity <- df2 %>%
  merge(never_smoker, by = "Gel_ID") %>%
  select(Smoking.status, DBS2.x) %>%
  rename(DBS2 = DBS2.x)

DBS2_activity_smoking_status <- data.frame(rbind(current_smoker_DBS2_activity, ex_smoker_DBS2_activity, never_smoker_DBS2_activity))
colnames(DBS2_activity_smoking_status) <- c("Smoking_status", "DBS2")
DBS2_activity_smoking_status$Smoking_status <- as.factor(DBS2_activity_smoking_status$Smoking_status)
DBS2_activity_smoking_status$DBS2 <- as.numeric(DBS2_activity_smoking_status$DBS2)
DBS2_activity_smoking_status$Smoking_status <- factor(DBS2_activity_smoking_status$Smoking_status, levels = c("Current", "Ex", "Never"))

my_comparisons <- split(t(combn(levels(DBS2_activity_smoking_status$Smoking_status), 2)), c(1,3,2)) #all pairwise combinations manual order

DBS2_smoking_status_stat <- pairwise.wilcox.test(DBS2_activity_smoking_status$DBS2, DBS2_activity_smoking_status$Smoking_status, p.adjust.method = "BH")

DBS2_smoking_status_boxplot <- ggplot(DBS2_activity_smoking_status, aes(x=Smoking_status, y=DBS2)) + 
  geom_boxplot(outlier.shape = NA, fill = "lightgrey") +
  geom_jitter(aes(colour = Smoking_status), size = 3) +
  theme_classic() +
  labs(x = "Smoking status", y="DBS2 mutations") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 20), 
        legend.title = element_text(face = "bold", size = 15),legend.text = element_text(size = 15),
        axis.title.x = element_text(face = "bold", size = 20), axis.text.x = element_text(size = 15), 
        axis.title.y = element_text(face = "bold", size = 20), axis.text.y = element_text(size = 15)) +
  scale_y_continuous(limits = c(0,6000), breaks = seq(0, 6000, by=1000)) +
  scale_colour_manual(values = c("Current" = "firebrick", "Ex" = "steelblue", 
                                 "Never" = "forestgreen"), name = "Smoking status") +
  geom_signif(comparisons = my_comparisons, textsize=6, fontface = "bold", step_increase = 0.1, map_signif_level = F,
              annotations = c(paste("p =", paste(signif(DBS2_smoking_status_stat$p.value[1], digits = 2))),
                              paste("p =", paste(signif(DBS2_smoking_status_stat$p.value[4], digits = 2))),
                              paste("p =", paste(signif(DBS2_smoking_status_stat$p.value[2], digits = 2))))) 

ggsave2(filename = "Figure_3E.pdf", plot = DBS2_smoking_status_boxplot, height = 8, width = 8, units = "in", path =  "~/file/path/Figures_final/", device = cairo_pdf)

# ID3
df3 <- COSMIC_IDS %>%
  tibble::rownames_to_column(var="Gel_ID") 
df3$Gel_ID <- as.numeric(df3$Gel_ID)

current_smoker_ID3_activity <- df3 %>%
  merge(current_smoker, by = "Gel_ID") %>%
  select(Smoking.status, ID3.x) %>%
  rename(ID3 = ID3.x)

ex_smoker_ID3_activity <- df3 %>%
  merge(ex_smoker, by = "Gel_ID") %>%
  select(Smoking.status, ID3.x) %>%
  rename(ID3 = ID3.x)

never_smoker_ID3_activity <- df3 %>%
  merge(never_smoker, by = "Gel_ID") %>%
  select(Smoking.status, ID3.x) %>%
  rename(ID3 = ID3.x)

ID3_activity_smoking_status <- data.frame(rbind(current_smoker_ID3_activity, ex_smoker_ID3_activity, never_smoker_ID3_activity))
colnames(ID3_activity_smoking_status) <- c("Smoking_status", "ID3")
ID3_activity_smoking_status$Smoking_status <- as.factor(ID3_activity_smoking_status$Smoking_status)
ID3_activity_smoking_status$ID3 <- as.numeric(ID3_activity_smoking_status$ID3)
ID3_activity_smoking_status$Smoking_status <- factor(ID3_activity_smoking_status$Smoking_status, levels = c("Current", "Ex", "Never"))

my_comparisons <- split(t(combn(levels(ID3_activity_smoking_status$Smoking_status), 2)), c(1,3,2)) #all pairwise combinations manual order

ID3_smoking_status_stat <- pairwise.wilcox.test(ID3_activity_smoking_status$ID3, ID3_activity_smoking_status$Smoking_status, p.adjust.method = "BH")

ID3_smoking_status_boxplot <- ggplot(ID3_activity_smoking_status, aes(x=Smoking_status, y=ID3)) + 
  geom_boxplot(outlier.shape = NA, fill = "lightgrey") +
  geom_jitter(aes(colour = Smoking_status), size = 3) +
  theme_classic() +
  labs(x = "Smoking status", y="ID3 mutations") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 20), 
        legend.title = element_text(face = "bold", size = 15),legend.text = element_text(size = 15),
        axis.title.x = element_text(face = "bold", size = 20), axis.text.x = element_text(size = 15), 
        axis.title.y = element_text(face = "bold", size = 20), axis.text.y = element_text(size = 15)) +
  scale_y_continuous(limits = c(0,11000), breaks = seq(0, 11000, by=1000)) +
  scale_colour_manual(values = c("Current" = "firebrick", "Ex" = "steelblue", 
                                 "Never" = "forestgreen"), name = "Smoking status") +
  geom_signif(comparisons = my_comparisons, textsize=6, fontface = "bold", step_increase = 0.1, map_signif_level = F,
              annotations = c(paste("p =", paste(signif(ID3_smoking_status_stat$p.value[1], digits = 2))),
                              paste("p =", paste(signif(ID3_smoking_status_stat$p.value[4], digits = 2))),
                              paste("p =", paste(signif(ID3_smoking_status_stat$p.value[2], digits = 2))))) 

ggsave2(filename = "Figure_3F.pdf", plot = ID3_smoking_status_boxplot, height = 8, width = 8, units = "in", path =  "~/file/path/Figures_final/", device = cairo_pdf)

# Figure 4 -------------------------------------------------------------------------
# Stats
age_at_cessation <- ex_smoker$Age.at.diagnosis - (ex_smoker$Time.as.ex.smoker..months./12)
summary(age_at_cessation)

duration_cessation <- ex_smoker$Time.as.ex.smoker..months./12
summary(duration_cessation)

# Figure 4A. Smoking cessation and SBS4 mutations
COSMIC_SBS_new <- COSMIC_SBS
colnames(COSMIC_SBS_new) <- paste0(colnames(COSMIC_SBS_new), "_counts")
lung_patients_combined_SBS_counts <- cbind(lung_patients_combined, COSMIC_SBS_new)

ever_smoker_SBS_counts <- lung_patients_combined_SBS_counts[,c("Time.as.ex.smoker..months.","SBS4_counts")]
colnames(ever_smoker_SBS_counts) <- c("Time", "SBS4_counts")

ever_smoker_SBS_counts_bins <- cut(x = ever_smoker_SBS_counts$Time, breaks = c(-Inf,0,120,240,360,480,600), include.lowest = T)
ever_smoker_SBS_counts_bins_combined <- cbind(ever_smoker_SBS_counts, ever_smoker_SBS_counts_bins)
ever_smoker_SBS_counts_bins_combined <- ever_smoker_SBS_counts_bins_combined[!is.na(ever_smoker_SBS_counts_bins_combined$Time),]

ever_smoker_SBS_counts_bins_combined_stat <- summary(lm(ever_smoker_SBS_counts_bins_combined$SBS4_counts ~ ever_smoker_SBS_counts_bins_combined$Time))

ever_smoker_SBS_counts_bins_combined_boxplot <- ggplot(data = ever_smoker_SBS_counts_bins_combined, aes(x = ever_smoker_SBS_counts_bins, y = SBS4_counts)) +
  geom_boxplot(colour = "black",  outlier.shape = NA, alpha = 0.75, fill = "firebrick") +
  geom_jitter(height = 0, width = 0.25) +
  theme_classic() +
  labs(x = "Duration of smoking cessation (years)", y="SBS4 mutations") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 20), legend.position = "none") +
  theme(axis.title.x = element_text(face = "bold", size = 20)) +
  theme(axis.text.x = element_text(size = 15)) +
  theme(axis.title.y = element_text(face = "bold", size = 20)) +
  theme(axis.text.y = element_text(size = 15)) +
  scale_y_continuous(limits = c(0,160000), breaks = seq(0, 160000, by=20000)) +
  scale_x_discrete(labels = c("Current smoker", "<10", "10-20", "20-30", "30-40", "40-50")) +  
  annotate("text", x = 5.5, y = 155000, size = 6, label = paste("R\u00b2 adj =", paste(signif(ever_smoker_SBS_counts_bins_combined_stat$adj.r.squared, digits = 2)))) +
  annotate("text", x = 5.5, y = 147500, size = 6, label = paste("p =", paste(signif(ever_smoker_SBS_counts_bins_combined_stat$coefficients[2,4], digits = 2))))

ggsave2(filename = "Figure_4A.pdf", plot = ever_smoker_SBS_counts_bins_combined_boxplot, height = 8, width = 8, units = "in", path = "~/file/path/Figures_final/", device = cairo_pdf)

# Boxplot with ANOVA
anova_SBS4 <- aov(SBS4_counts ~ ever_smoker_SBS_counts_bins, data = ever_smoker_SBS_counts_bins_combined)
summary(anova_SBS4)
anova_SBS4_p_value <- summary(anova_SBS4)[[1]]$`Pr(>F)`[1]
formatted_p_value <- paste("p=", signif(anova_SBS4_p_value, digits=3))

ever_smoker_SBS_counts_bins_combined_boxplot_ANOVA <- ggplot(data = ever_smoker_SBS_counts_bins_combined, aes(x = ever_smoker_SBS_counts_bins, y = SBS4_counts)) +
  geom_boxplot(colour = "black",  outlier.shape = NA, alpha = 0.75, fill = "firebrick") +
  geom_jitter(height = 0, width = 0.25) +
  theme_classic() +
  labs(x = "Duration of smoking cessation (years)", y="SBS4 mutations") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 20), legend.position = "none") +
  theme(axis.title.x = element_text(face = "bold", size = 20)) +
  theme(axis.text.x = element_text(size = 15)) +
  theme(axis.title.y = element_text(face = "bold", size = 20)) +
  theme(axis.text.y = element_text(size = 15)) +
  scale_y_continuous(limits = c(0,160000), breaks = seq(0, 160000, by=20000)) +
  scale_x_discrete(labels = c("Current smoker", "<10", "10-20", "20-30", "30-40", "40-50")) + 
  annotate("text", x = 5.5, y = 147500, size = 6, label = formatted_p_value)

ggsave2(filename = "Figure_4A3.pdf", plot = ever_smoker_SBS_counts_bins_combined_boxplot_ANOVA, height = 8, width = 8, units = "in", path = "~/file/path/Figures_final/", device = cairo_pdf)

# do with linear regression instead of box plots 
COSMIC_SBS_new <- COSMIC_SBS
colnames(COSMIC_SBS_new) <- paste0(colnames(COSMIC_SBS_new), "_counts")
lung_patients_combined_SBS_counts <- cbind(lung_patients_combined, COSMIC_SBS_new)

ever_smoker_SBS_counts <- lung_patients_combined_SBS_counts[,c("Time.as.ex.smoker..months.","SBS4_counts")]
colnames(ever_smoker_SBS_counts) <- c("Time", "SBS4_counts")
# SBS4 
ever_smoker_SBS_counts_bins_combined_linear <- ggplot(data = ever_smoker_SBS_counts_bins_combined, aes(x = Time, y = SBS4_counts)) +
  geom_point()+
  geom_smooth(method="lm", se=TRUE, color="steelblue")+
  theme_classic() +
  labs(x = "Duration of smoking cessation (months)", y="SBS4 mutations") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 20), legend.position = "none") +
  theme(axis.title.x = element_text(face = "bold", size = 20)) +
  theme(axis.text.x = element_text(size = 15)) +
  theme(axis.title.y = element_text(face = "bold", size = 20)) +
  theme(axis.text.y = element_text(size = 15)) + 
  scale_y_continuous(limits = c(0,160000), breaks = seq(0, 160000, by=20000)) + 
  annotate("text", x = 100, y = 155000, size = 6, label = paste("R\u00b2 adj =", paste(signif(ever_smoker_SBS_counts_bins_combined_stat$adj.r.squared, digits = 2)))) +
  annotate("text", x = 100, y = 147500, size = 6, label = paste("p =", paste(signif(ever_smoker_SBS_counts_bins_combined_stat$coefficients[2,4], digits = 2))))

ggsave2(filename = "Figure_4A2.pdf", plot = ever_smoker_SBS_counts_bins_combined_linear, height = 8, width = 8, units = "in", path = "~/file/path/Figures_final/", device = cairo_pdf)

#Figure 4B
COSMIC_DBS_new <- COSMIC_DBS
colnames(COSMIC_DBS_new) <- paste0(colnames(COSMIC_DBS_new), "_counts")
lung_patients_combined_DBS_counts <- cbind(lung_patients_combined, COSMIC_DBS_new)

ever_smoker_DBS_counts <- lung_patients_combined_DBS_counts[,c("Time.as.ex.smoker..months.","DBS2_counts")]
colnames(ever_smoker_DBS_counts) <- c("Time", "DBS2_counts")

ever_smoker_DBS_counts_bins <- cut(x = ever_smoker_DBS_counts$Time, breaks = c(-Inf,0,120,240,360,480,600), include.lowest = T)
ever_smoker_DBS_counts_bins_combined <- cbind(ever_smoker_DBS_counts, ever_smoker_DBS_counts_bins)
ever_smoker_DBS_counts_bins_combined <- ever_smoker_DBS_counts_bins_combined[!is.na(ever_smoker_DBS_counts_bins_combined$Time),]

ever_smoker_DBS_counts_bins_combined_stat <- summary(lm(ever_smoker_DBS_counts_bins_combined$DBS2_counts ~ ever_smoker_DBS_counts_bins_combined$Time))

ever_smoker_DBS_counts_bins_combined_boxplot <- ggplot(data = ever_smoker_DBS_counts_bins_combined, aes(x = ever_smoker_DBS_counts_bins, y = DBS2_counts)) +
  geom_boxplot(colour = "black",  outlier.shape = NA, alpha = 0.75, fill = "darkred") +
  geom_jitter(height = 0, width = 0.25) +
  theme_classic() +
  labs(x = "Duration of smoking cessation (years)", y="DBS2 mutations") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 20), legend.position = "none") +
  theme(axis.title.x = element_text(face = "bold", size = 20)) +
  theme(axis.text.x = element_text(size = 15)) +
  theme(axis.title.y = element_text(face = "bold", size = 20)) +
  theme(axis.text.y = element_text(size = 15)) +
  scale_y_continuous(limits = c(0,4500), breaks = seq(0, 4500, by=500)) +
  scale_x_discrete(labels = c("Current smoker", "<10", "10-20", "20-30", "30-40", "40-50")) +  
  annotate("text", x = 5.5, y = 4300, size = 6, label = paste("R\u00b2 adj =", paste(signif(ever_smoker_DBS_counts_bins_combined_stat$adj.r.squared, digits = 2)))) +
  annotate("text", x = 5.5, y = 4000, size = 6, label = paste("p =", paste(signif(ever_smoker_DBS_counts_bins_combined_stat$coefficients[2,4], digits = 2))))

ggsave2(filename = "Figure_4B.pdf", plot = ever_smoker_DBS_counts_bins_combined_boxplot, height = 8, width = 8, units = "in", path = "~/file/path/Figures_final/", device = cairo_pdf)

# Box plot with ANOVA 
anova_DBS2 <- aov(DBS2_counts ~ ever_smoker_DBS_counts_bins, data = ever_smoker_DBS_counts_bins_combined)
summary(anova_DBS2)
anova_DBS2_p_value <- summary(anova_DBS2)[[1]]$`Pr(>F)`[1]
formatted_p_value <- paste("p=", signif(anova_DBS2_p_value, digits=3))

ever_smoker_DBS_counts_bins_combined_boxplot_ANOVA <- ggplot(data = ever_smoker_DBS_counts_bins_combined, aes(x = ever_smoker_DBS_counts_bins, y = DBS2_counts)) +
  geom_boxplot(colour = "black",  outlier.shape = NA, alpha = 0.75, fill = "darkred") +
  geom_jitter(height = 0, width = 0.25) +
  theme_classic() +
  labs(x = "Duration of smoking cessation (years)", y="DBS2 mutations") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 20), legend.position = "none") +
  theme(axis.title.x = element_text(face = "bold", size = 20)) +
  theme(axis.text.x = element_text(size = 15)) +
  theme(axis.title.y = element_text(face = "bold", size = 20)) +
  theme(axis.text.y = element_text(size = 15)) +
  scale_y_continuous(limits = c(0,4500), breaks = seq(0, 4500, by=500)) +
  scale_x_discrete(labels = c("Current smoker", "<10", "10-20", "20-30", "30-40", "40-50")) +  
  annotate("text", x = 5.5, y = 4000, size = 6, label = formatted_p_value)

ggsave2(filename = "Figure_4B3.pdf", plot = ever_smoker_DBS_counts_bins_combined_boxplot_ANOVA, height = 8, width = 8, units = "in", path = "~/file/path/Figures_final/", device = cairo_pdf)

# do with linear regression instead of box plots 
# DBS2
ever_smoker_DBS_counts_bins_combined_linear <- ggplot(data = ever_smoker_DBS_counts_bins_combined, aes(x = Time, y = DBS2_counts)) +
  geom_point()+
  geom_smooth(method="lm", se=TRUE, color="steelblue")+
  theme_classic() +
  labs(x = "Duration of smoking cessation (months)", y="DBS2 mutations") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 20), legend.position = "none") +
  theme(axis.title.x = element_text(face = "bold", size = 20)) +
  theme(axis.text.x = element_text(size = 15)) +
  theme(axis.title.y = element_text(face = "bold", size = 20)) +
  theme(axis.text.y = element_text(size = 15)) +
  scale_y_continuous(limits = c(0,4500), breaks = seq(0, 4500, by=500)) +
  annotate("text", x = 100, y = 4300, size = 6, label = paste("R\u00b2 adj =", paste(signif(ever_smoker_DBS_counts_bins_combined_stat$adj.r.squared, digits = 2)))) +
  annotate("text", x = 100, y = 4000, size = 6, label = paste("p =", paste(signif(ever_smoker_DBS_counts_bins_combined_stat$coefficients[2,4], digits = 2))))

ggsave2(filename = "Figure_4B2.pdf", plot = ever_smoker_DBS_counts_bins_combined_linear, height = 8, width = 8, units = "in", path = "~/file/path/Figures_final/", device = cairo_pdf)

#Figure 4C
COSMIC_IDS_new <- COSMIC_IDS
colnames(COSMIC_IDS_new) <- paste0(colnames(COSMIC_IDS_new), "_counts")
lung_patients_combined_IDS_counts <- cbind(lung_patients_combined, COSMIC_IDS_new)

ever_smoker_IDS_counts <- lung_patients_combined_IDS_counts[,c("Time.as.ex.smoker..months.","ID3_counts")]
colnames(ever_smoker_IDS_counts) <- c("Time", "ID3_counts")

ever_smoker_IDS_counts_bins <- cut(x = ever_smoker_DBS_counts$Time, breaks = c(-Inf,0,120,240,360,480,600), include.lowest = T)
ever_smoker_IDS_counts_bins_combined <- cbind(ever_smoker_IDS_counts, ever_smoker_IDS_counts_bins)
ever_smoker_IDS_counts_bins_combined <- ever_smoker_IDS_counts_bins_combined[!is.na(ever_smoker_IDS_counts_bins_combined$Time),]

ever_smoker_IDS_counts_bins_combined_stat <- summary(lm(ever_smoker_IDS_counts_bins_combined$ID3_counts ~ ever_smoker_IDS_counts_bins_combined$Time))

ever_smoker_IDS_counts_bins_combined_boxplot <- ggplot(data = ever_smoker_IDS_counts_bins_combined, aes(x = ever_smoker_IDS_counts_bins, y = ID3_counts)) +
  geom_boxplot(colour = "black",  outlier.shape = NA, alpha = 0.75, fill = "bisque2") +
  geom_jitter(height = 0, width = 0.25) +
  theme_classic() +
  labs(x = "Duration of smoking cessation (years)", y="ID3 mutations") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 20), legend.position = "none") +
  theme(axis.title.x = element_text(face = "bold", size = 20)) +
  theme(axis.text.x = element_text(size = 15)) +
  theme(axis.title.y = element_text(face = "bold", size = 20)) +
  theme(axis.text.y = element_text(size = 15)) +
  scale_y_continuous(limits = c(0,8000), breaks = seq(0, 8000, by=1000)) +
  scale_x_discrete(labels = c("Current smoker", "<10", "10-20", "20-30", "30-40", "40-50")) +  
  annotate("text", x = 5.5, y = 7800, size = 6, label = paste("R\u00b2 adj =", paste(signif(ever_smoker_IDS_counts_bins_combined_stat$adj.r.squared, digits = 2)))) +
  annotate("text", x = 5.5, y = 7500, size = 6, label = paste("p =", paste(signif(ever_smoker_IDS_counts_bins_combined_stat$coefficients[2,4], digits = 2))))

ggsave2(filename = "Figure_4C.pdf", plot = ever_smoker_IDS_counts_bins_combined_boxplot, height = 8, width = 8, units = "in", path = "~/file/path/Figures_final/", device = cairo_pdf)

# plot with anova
anova_ID3 <- aov(ID3_counts ~ ever_smoker_IDS_counts_bins, data = ever_smoker_IDS_counts_bins_combined)
summary(anova_ID3)
anova_ID3_p_value <- summary(anova_ID3)[[1]]$`Pr(>F)`[1]
formatted_p_value <- paste("p=", signif(anova_ID3_p_value, digits=3))

ever_smoker_IDS_counts_bins_combined_boxplot_ANOVA <- ggplot(data = ever_smoker_IDS_counts_bins_combined, aes(x = ever_smoker_IDS_counts_bins, y = ID3_counts)) +
  geom_boxplot(colour = "black",  outlier.shape = NA, alpha = 0.75, fill = "bisque2") +
  geom_jitter(height = 0, width = 0.25) +
  theme_classic() +
  labs(x = "Duration of smoking cessation (years)", y="ID3 mutations") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 20), legend.position = "none") +
  theme(axis.title.x = element_text(face = "bold", size = 20)) +
  theme(axis.text.x = element_text(size = 15)) +
  theme(axis.title.y = element_text(face = "bold", size = 20)) +
  theme(axis.text.y = element_text(size = 15)) +
  scale_y_continuous(limits = c(0,8000), breaks = seq(0, 8000, by=1000)) +
  scale_x_discrete(labels = c("Current smoker", "<10", "10-20", "20-30", "30-40", "40-50")) +  
  annotate("text", x = 5.5, y = 7500, size = 6, label = formatted_p_value)

ggsave2(filename = "Figure_4C3.pdf", plot = ever_smoker_IDS_counts_bins_combined_boxplot_ANOVA, height = 8, width = 8, units = "in", path = "~/file/path/Figures_final/", device = cairo_pdf)

# do with linear regression instead of box plots 
# ID3 
ever_smoker_IDS_counts_bins_combined_linear <- ggplot(data = ever_smoker_IDS_counts_bins_combined, aes(x = Time, y = ID3_counts)) +
  geom_point()+
  geom_smooth(method="lm", se=TRUE, color="steelblue")+
  theme_classic() +
  labs(x = "Duration of smoking cessation (months)", y="ID3 mutations") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 20), legend.position = "none") +
  theme(axis.title.x = element_text(face = "bold", size = 20)) +
  theme(axis.text.x = element_text(size = 15)) +
  theme(axis.title.y = element_text(face = "bold", size = 20)) +
  theme(axis.text.y = element_text(size = 15)) +
  scale_y_continuous(limits = c(0,8000), breaks = seq(0, 8000, by=1000)) +
  annotate("text", x = 100, y = 7800, size = 6, label = paste("R\u00b2 adj =", paste(signif(ever_smoker_IDS_counts_bins_combined_stat$adj.r.squared, digits = 2)))) +
  annotate("text", x = 100, y = 7500, size = 6, label = paste("p =", paste(signif(ever_smoker_IDS_counts_bins_combined_stat$coefficients[2,4], digits = 2))))

ggsave2(filename = "Figure_4C2.pdf", plot = ever_smoker_IDS_counts_bins_combined_linear, height = 8, width = 8, units = "in", path = "~/file/path/Figures_final/", device = cairo_pdf)

# Cessation plot combined 
cessation_combined <- plot_grid(ever_smoker_SBS_counts_bins_combined_boxplot, ever_smoker_DBS_counts_bins_combined_boxplot, ever_smoker_IDS_counts_bins_combined_boxplot, Cessation_time_TMB_bins_combined_boxplot,align = "h", nrow = 1)

ggsave2(filename = "Figure_4AtoD.pdf", plot = cessation_combined, height = 6, width = 26, units = "in", path = "~/file/path/Figures_final/", device = cairo_pdf)

# Figure 4D - smoking cessation, TMB
Cessation_time_TMB <- lung_patients_combined[,c("Time.as.ex.smoker..months.","TMB")]
colnames(Cessation_time_TMB) <- c("Time", "TMB")
#ex_smoker_time_TMB$Time <- (ex_smoker_time_TMB$Time)/12 convert to years

Cessation_time_TMB <- Cessation_time_TMB[!is.na(Cessation_time_TMB$Time),]
Cessation_time_TMB_bins <- cut(x = Cessation_time_TMB$Time, breaks = c(-Inf,0,120,240,360,480,600), include.lowest = T)
Cessation_time_TMB_bins_combined <- cbind(Cessation_time_TMB, Cessation_time_TMB_bins)

Cessation_time_TMB_bins_combined$Cessation_time_TMB_bins <-as.factor(Cessation_time_TMB_bins_combined$Cessation_time_TMB_bins)
Cessation_time_TMB_bins_combined$Cessation_time_TMB_bins <- factor(Cessation_time_TMB_bins_combined$Cessation_time_TMB_bins, 
                                                                   levels = c("[-Inf,0]", "(0,120]", "(120,240]", "(240,360]", "(360,480]", "(480,600]"))
Cessation_time_TMB_bins_combined_stat <- summary(lm(Cessation_time_TMB_bins_combined$TMB ~ Cessation_time_TMB_bins_combined$Time))

Cessation_time_TMB_bins_combined_boxplot <- ggplot(data = Cessation_time_TMB_bins_combined, aes(x = Cessation_time_TMB_bins, y = TMB)) +
  geom_boxplot(colour = "black",  outlier.shape = NA, fill = "lightgrey") +
  geom_jitter(height = 0, width = 0.25) +
  theme_classic() +
  labs(x = "Duration of smoking cessation (years)", y="TMB (mutations/Mb)") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 20), legend.position = "none") +
  theme(axis.title.x = element_text(face = "bold", size = 20)) +
  theme(axis.text.x = element_text(size = 15)) +
  theme(axis.title.y = element_text(face = "bold", size = 20)) +
  theme(axis.text.y = element_text(size = 15)) +
  scale_y_continuous(limits = c(0,50), breaks = seq(0, 50, by=5)) +
  scale_x_discrete(labels = c("Current smoker", "<10", "10-20", "20-30", "30-40", "40-50")) +  
  annotate("text", x = 5.5, y = 50, size = 6, label = paste("R\u00b2 adj =", paste(signif(Cessation_time_TMB_bins_combined_stat$adj.r.squared, digits = 2)))) +
  annotate("text", x = 5.5, y = 48, size = 6, label = paste("p =", paste(signif(Cessation_time_TMB_bins_combined_stat$coefficients[2,4], digits = 2))))

ggsave2(filename = "Figure_4D.pdf", plot = Cessation_time_TMB_bins_combined_boxplot, height = 8, width = 8, units = "in", path = "~/file/path/Figures_final/", device = cairo_pdf)

# plot with anova

anova_TMB <- aov(TMB~ Cessation_time_TMB_bins, data = Cessation_time_TMB_bins_combined)
summary(anova_TMB)
anova_TMB_p_value <- summary(anova_TMB)[[1]]$`Pr(>F)`[1]
formatted_p_value <- paste("p=", signif(anova_TMB_p_value, digits=3))

Cessation_time_TMB_bins_combined_boxplot_ANOVA <- ggplot(data = Cessation_time_TMB_bins_combined, aes(x = Cessation_time_TMB_bins, y = TMB)) +
  geom_boxplot(colour = "black",  outlier.shape = NA, fill = "lightgrey") +
  geom_jitter(height = 0, width = 0.25) +
  theme_classic() +
  labs(x = "Duration of smoking cessation (years)", y="TMB (mutations/Mb)") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 20), legend.position = "none") +
  theme(axis.title.x = element_text(face = "bold", size = 20)) +
  theme(axis.text.x = element_text(size = 15)) +
  theme(axis.title.y = element_text(face = "bold", size = 20)) +
  theme(axis.text.y = element_text(size = 15)) +
  scale_y_continuous(limits = c(0,50), breaks = seq(0, 50, by=5)) +
  scale_x_discrete(labels = c("Current smoker", "<10", "10-20", "20-30", "30-40", "40-50")) +  
  annotate("text", x = 5.5, y = 48, size = 6, label = formatted_p_value)

ggsave2(filename = "Figure_4D3.pdf", plot = Cessation_time_TMB_bins_combined_boxplot_ANOVA, height = 8, width = 8, units = "in", path = "~/file/path/Figures_final/", device = cairo_pdf)


# TMB linear
Cessation_time_TMB_bins_combined_linear <- ggplot(data = Cessation_time_TMB_bins_combined, aes(x = Time, y = TMB)) +
  geom_point()+
  geom_smooth(method="lm", se=TRUE, color="steelblue")+
  theme_classic() +
  labs(x = "Duration of smoking cessation (months)", y="TMB (mutations/Mb)") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 20), legend.position = "none") +
  theme(axis.title.x = element_text(face = "bold", size = 20)) +
  theme(axis.text.x = element_text(size = 15)) +
  theme(axis.title.y = element_text(face = "bold", size = 20)) +
  theme(axis.text.y = element_text(size = 15)) +
  scale_y_continuous(limits = c(0,50), breaks = seq(0, 50, by=5)) +
  annotate("text", x = 100, y = 50, size = 6, label = paste("R\u00b2 adj =", paste(signif(Cessation_time_TMB_bins_combined_stat$adj.r.squared, digits = 2)))) +
  annotate("text", x = 100, y = 48, size = 6, label = paste("p =", paste(signif(Cessation_time_TMB_bins_combined_stat$coefficients[2,4], digits = 2))))

ggsave2(filename = "Figure_4D2.pdf", plot = Cessation_time_TMB_bins_combined_linear, height = 8, width = 8, units = "in", path = "~/file/path/Figures_final/", device = cairo_pdf)


# Figure 4E - ex smoker venn diagram
venn_df_ex <- venn_df %>% 
  filter(Smoking.status == "Ex")

area1 <- venn_df_ex %>%
  filter(!is.na(SBS4) & is.na(DBS2) & is.na(ID3)) %>%
  nrow()
area2 <- venn_df_ex %>%
  filter(is.na(SBS4) & !is.na(DBS2) & is.na(ID3)) %>%
  nrow()
area3 <- venn_df_ex %>%
  filter(is.na(SBS4) & is.na(DBS2) & !is.na(ID3)) %>%
  nrow()

n12 <-nrow(unique(venn_df_ex[which(!is.na(venn_df_ex$SBS4) & !is.na(venn_df_ex$DBS2) 
                                   & is.na(venn_df_ex$ID3)),]))
n13 <- nrow(unique(venn_df_ex[which(!is.na(venn_df_ex$SBS4) & is.na(venn_df_ex$DBS2) 
                                    & !is.na(venn_df_ex$ID3)),]))
n23 <- nrow(unique(venn_df_ex[which(is.na(venn_df_ex$SBS4) & !is.na(venn_df_ex$DBS2) 
                                    & !is.na(venn_df_ex$ID3)),]))

n123 <- nrow(unique(venn_df_ex[which(!is.na(venn_df_ex$SBS4) & !is.na(venn_df_ex$DBS2) 
                                     & !is.na(venn_df_ex$ID3)),]))

fit4 <- c("area1" = area1, "area2" = area2, "area3" = area3, "area1&area2" = n12, "area1&area3" = n13, 
          "area2&area3" = n23, "area1&area2&area3" = n123)

venn4 <- venn(fit4)
eulerr_options(labels = list(cex = 2))
euler_4 <- plot(venn4, labels = c("SBS4", "DBS2", "ID3"), quantities = list(cex = 2), fills = c("grey", "peachpuff", "lightblue"))

ggsave2(filename = "Figure_4E.pdf", plot = euler_4, height = 8, width = 8, units = "in", path = "~/file/path/Figures_final/", device = cairo_pdf)

# Figure 4F
#Survival plots in ever-smokers using quartiles
#Pack-years in ex-smokers according to time of cessation
pack_year_cessation <- lung_patients[,c("Smoking..pack.years.", "Time.as.ex.smoker..months.")]
colnames(pack_year_cessation) <- c("Pack_years", "Time")
pack_year_cessation <- pack_year_cessation[!is.na(pack_year_cessation$Time),]
summary(pack_year_cessation$Pack_years)

Pack_year_cessation_bins <- cut(x = pack_year_cessation$Pack_years, breaks = c(2.5,20,35,50,165), include.lowest = T)
Pack_year_cessation_bins_combined <- cbind(pack_year_cessation, Pack_year_cessation_bins)
Pack_year_cessation_bins_combined$status <- 1
colnames(Pack_year_cessation_bins_combined) <- c("Pack_years", "Time", "Pack_year_bins", "status")
table(Pack_year_cessation_bins_combined$Pack_year_bins, useNA = "always")

res.cox_4 <- coxph(Surv(Time, status) ~ Pack_year_bins, data = Pack_year_cessation_bins_combined)
summary(res.cox_4)
fit <- survfit(Surv(Time, status) ~ Pack_year_bins, data = Pack_year_cessation_bins_combined)
surv_summary(fit)
log_rank_p <- surv_pvalue(fit)

log_rank_p_pairwise <- pairwise_survdiff(Surv(Time, status) ~ Pack_year_bins, data = Pack_year_cessation_bins_combined, p.adjust.method = "BH")

Cumulative_incidence_pack_year_plot <- ggsurvplot(fit, data = Pack_year_cessation_bins_combined, fun = "event", 
                                                  ggtheme = theme_classic(), xlab = "Smoking cessation time (months)", ylab = "Cumulative NSCLC incidence", 
                                                  palette = "jco", combine = T, conf.int = T, 
                                                  pval = F, 
                                                  surv.median.line = "hv", risk.table = F, cumevents = T, tables.height = 0.3, cumevents.title = "",
                                                  legend.title = "Pack-years", legend.labs = c("2.5-20", "21-35", "36-50", "51-165"))

Cumulative_incidence_pack_year_plot$plot <- Cumulative_incidence_pack_year_plot$plot +
  theme(legend.title = element_text(size = 20, face = "bold"), legend.text = element_text(size = 20), 
        axis.text = element_text(size = 20), axis.title = element_text(size = 20, face = "bold")) +
  scale_y_continuous(limits = c(0,1), breaks = seq(0,1,by = 0.25), expand = c(0,0)) +
  scale_x_continuous(limits = c(0,600), breaks = seq(0,600,by = 100), expand = c(0,15)) +
  annotate("text", size = 6, x = 500, y = 0.05, label = paste("p =", paste(signif(log_rank_p$pval, digits = 2)), "by log-rank test"))

Cumulative_incidence_pack_year_plot$cumevents <- Cumulative_incidence_pack_year_plot$cumevents +
  theme(axis.text = element_text(size = 20), axis.title = element_text(size = 20, face = "bold")) +
  scale_x_continuous(limits = c(0,600), breaks = seq(0,600,by = 100), expand = c(0,15))

Cumulative_incidence_pack_year_plot$cumevents$layers[[1]]$aes_params$size <- 6

pdf(file = "~/file/path/Figures_final/Figure_4F.pdf",
    width = 20, height = 12)
print(Cumulative_incidence_pack_year_plot, newpage = FALSE)
dev.off()

# Supp Figure 3 #####
# Smoking abstinence before age 30
#Ex-smokers who quit under 30 years old
ex_smoker_ages <- ex_smoker[,c("Gel_ID", "Age.at.diagnosis", "TMB", "Smoking..pack.years.", "Time.as.ex.smoker..months.")]
ex_smoker_ages$cessation_age <- ex_smoker_ages$Age.at.diagnosis - (ex_smoker_ages$Time.as.ex.smoker..months./12)
colnames(ex_smoker_ages) <- c("Gel_ID", "Age", "TMB", "Py", "Cessation_duration", "Cessation_age")
ex_smoker_ages_under_30 <- data.frame(ex_smoker_ages[which(ex_smoker_ages$Cessation_age < 30),])
summary(ex_smoker_ages_under_30$Py)
summary(ex_smoker_ages_under_30$TMB)

# Ex smokers over 30
ex_smoker_ages_over_30 <- data.frame(ex_smoker_ages[which(ex_smoker_ages$Cessation_age >= 30),])
nrow(ex_smoker_ages_over_30)
summary(ex_smoker_ages_over_30$TMB)
summary(never_smoker$TMB)

# Create a dataframe with all 3 groups (never, ex <30, ex >30)
TMB_ex30_never_smokers <- data.frame(rbind(cbind("Cessation under 30", ex_smoker_ages_under_30$TMB),
                                           cbind("Cessation over 30", ex_smoker_ages_over_30$TMB),
                                           cbind("Never smoker", never_smoker$TMB)))
colnames(TMB_ex30_never_smokers) <- c("smoking_age_cohort", "TMB")
TMB_ex30_never_smokers$smoking_age_cohort <- as.factor(TMB_ex30_never_smokers$smoking_age_cohort)
TMB_ex30_never_smokers$TMB <- as.numeric(TMB_ex30_never_smokers$TMB)
TMB_ex30_never_smokers$smoking_age_cohort <- factor(TMB_ex30_never_smokers$smoking_age_cohort, levels = c("Cessation over 30", 
                                                                                                          "Cessation under 30",
                                                                                                          "Never smoker"))
my_comparisons <- split(t(combn(levels(TMB_ex30_never_smokers$smoking_age_cohort), 2)), c(1,3,2)) #all pairwise combinations manual order

TMB_cessation_age_stat <- pairwise.wilcox.test(TMB_ex30_never_smokers$TMB, TMB_ex30_never_smokers$smoking_age_cohort, p.adjust.method = "BH")

# Fig S3A
TMB_cessation_age_boxplot <- ggplot(TMB_ex30_never_smokers, aes(x=smoking_age_cohort, y=TMB)) + 
  geom_boxplot(outlier.shape = NA, fill = "lightgrey") +
  geom_jitter(aes(colour = smoking_age_cohort), size = 3) +
  theme_classic() +
  labs(x = "Age at smoking cessation", y="TMB (mutations/Mb)") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 20), 
        legend.title = element_text(face = "bold", size = 15),legend.text = element_text(size = 15),
        axis.title.x = element_text(face = "bold", size = 20), axis.text.x = element_text(size = 15), 
        axis.title.y = element_text(face = "bold", size = 20), axis.text.y = element_text(size = 15)) +
  scale_y_continuous(limits = c(0,60), breaks = seq(0, 60, by=5)) +
  scale_colour_manual(values = c("Cessation over 30" = "firebrick", "Cessation under 30" = "steelblue", 
                                 "Never smoker" = "forestgreen"), name = " ") +
  geom_signif(comparisons = my_comparisons, textsize=6, fontface = "bold", step_increase = 0.1, map_signif_level = F,
              annotations = c(paste("p =", paste(signif(TMB_cessation_age_stat$p.value[1], digits = 2))),
                              paste("p =", paste(signif(TMB_cessation_age_stat$p.value[4], digits = 2))),
                              paste("p =", paste(signif(TMB_cessation_age_stat$p.value[2], digits = 2))))) 

ggsave2(filename = "Figure_S3A.pdf", plot = TMB_cessation_age_boxplot, height = 8, width =10, units = "in", path =  "~/file/path/Figures_final/", device = cairo_pdf)

# Figure S3B
ever_smoker_SBS_stacked_counts_bins_combined <- ever_smoker_SBS_counts_bins_combined
ever_smoker_SBS_stacked_counts_bins_combined$SBS4_counts[ever_smoker_SBS_stacked_counts_bins_combined$SBS4_counts > 0] <- "SBS4-positive"
ever_smoker_SBS_stacked_counts_bins_combined$SBS4_counts[ever_smoker_SBS_stacked_counts_bins_combined$SBS4_counts == 0] <- "SBS4-negative"

proportion_ever_smokers_SBS4 <- data.frame(with(ever_smoker_SBS_stacked_counts_bins_combined, table(SBS4_counts, ever_smoker_SBS_counts_bins)))

proportion_ever_smokers_SBS4_barplot <- ggplot(data = proportion_ever_smokers_SBS4, aes(x = ever_smoker_SBS_counts_bins, y = Freq, fill = SBS4_counts, beside = T)) +
  geom_bar(stat = "identity", colour = "black") +
  theme_classic() +
  labs(x = "Duration of smoking cessation (years)", y="Number of patients", fill = "SBS4 status") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 20)) +
  theme(axis.title.x = element_text(face = "bold", size = 20)) +
  theme(axis.text.x = element_text(size = 15)) +
  theme(axis.title.y = element_text(face = "bold", size = 20)) +
  theme(axis.text.y = element_text(size = 15)) +
  theme(legend.position = c(0.85,0.95), legend.title = element_text(face = "bold", size = 15), legend.text = element_text(size = 15)) +
  scale_y_continuous(limits = c(0,40), breaks = seq(0, 40, by=5)) +
  scale_x_discrete(labels = c("Current smokers", "<10", "11-20", "21-30", "31-40", "41-50")) +
  scale_fill_manual(values = c("steelblue", "firebrick"))

ggsave2(filename = "Figure_S3B.pdf", plot = proportion_ever_smokers_SBS4_barplot, height = 8, width = 8, units = "in", path = "~/file/path/Figures_final/", device = cairo_pdf)

# Figure S3C
#Proportion of ex-smokers who develop DBS2 positive and negative cancers
ever_smoker_DBS_stacked_counts_bins_combined <- ever_smoker_DBS_counts_bins_combined
ever_smoker_DBS_stacked_counts_bins_combined$DBS2_counts[ever_smoker_DBS_stacked_counts_bins_combined$DBS2_counts > 0] <- "DBS2-positive"
ever_smoker_DBS_stacked_counts_bins_combined$DBS2_counts[ever_smoker_DBS_stacked_counts_bins_combined$DBS2_counts == 0] <- "DBS2-negative"

proportion_ever_smokers_DBS2 <- data.frame(with(ever_smoker_DBS_stacked_counts_bins_combined, table(DBS2_counts, ever_smoker_DBS_counts_bins)))

proportion_ever_smokers_DBS2_barplot <- ggplot(data = proportion_ever_smokers_DBS2, aes(x = ever_smoker_DBS_counts_bins, y = Freq, fill = DBS2_counts, beside = T)) +
  geom_bar(stat = "identity", colour = "black") +
  theme_classic() +
  labs(x = "Duration of smoking cessation (years)", y="Number of patients", fill = "DBS2 status") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 20)) +
  theme(axis.title.x = element_text(face = "bold", size = 20)) +
  theme(axis.text.x = element_text(size = 15)) +
  theme(axis.title.y = element_text(face = "bold", size = 20)) +
  theme(axis.text.y = element_text(size = 15)) +
  theme(legend.position = c(0.85,0.95), legend.title = element_text(face = "bold", size = 15), legend.text = element_text(size = 15)) +
  scale_y_continuous(limits = c(0,40), breaks = seq(0, 40, by=5)) +
  scale_x_discrete(labels = c("Current smokers", "<10", "11-20", "21-30", "31-40", "41-50")) +
  scale_fill_manual(values = c("steelblue", "darkred"))

ggsave2(filename = "Figure_S3C.pdf", plot = proportion_ever_smokers_DBS2_barplot, height = 8, width = 8, units = "in", path = "~/file/path/Figures_final/", device = cairo_pdf)

# Figure S3D
# Proportion of ex-smokers who develop ID3 positive and negative cancers
ever_smoker_IDS_stacked_counts_bins_combined <- ever_smoker_IDS_counts_bins_combined
ever_smoker_IDS_stacked_counts_bins_combined$ID3_counts[ever_smoker_IDS_stacked_counts_bins_combined$ID3_counts > 0] <- "ID3-positive"
ever_smoker_IDS_stacked_counts_bins_combined$ID3_counts[ever_smoker_IDS_stacked_counts_bins_combined$ID3_counts == 0] <- "ID3-negative"

proportion_ever_smokers_ID3 <- data.frame(with(ever_smoker_IDS_stacked_counts_bins_combined, table(ID3_counts, ever_smoker_IDS_counts_bins)))

proportion_ever_smokers_ID3_barplot <- ggplot(data = proportion_ever_smokers_ID3, aes(x = ever_smoker_IDS_counts_bins, y = Freq, fill = ID3_counts, beside = T)) +
  geom_bar(stat = "identity", colour = "black") +
  theme_classic() +
  labs(x = "Duration of smoking cessation (years)", y="Number of patients", fill = "ID3 status") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 20)) +
  theme(axis.title.x = element_text(face = "bold", size = 20)) +
  theme(axis.text.x = element_text(size = 15)) +
  theme(axis.title.y = element_text(face = "bold", size = 20)) +
  theme(axis.text.y = element_text(size = 15)) +
  theme(legend.position = c(0.85,0.95), legend.title = element_text(face = "bold", size = 15), legend.text = element_text(size = 15)) +
  scale_y_continuous(limits = c(0,40), breaks = seq(0, 40, by=5)) +
  scale_x_discrete(labels = c("Current smokers", "<10", "11-20", "21-30", "31-40", "41-50")) +
  scale_fill_manual(values = c("steelblue", "bisque2"))

ggsave2(filename = "Figure_S3D.pdf", plot = proportion_ever_smokers_ID3_barplot, height = 8, width = 8, units = "in", path = "~/file/path/Figures_final/", device = cairo_pdf)

# Figure S3E 
# Produce table with median CFI 
median_time_per_group <- Pack_year_cessation_bins_combined %>%
  group_by(Pack_year_bins) %>%
  summarise(median_time = median(Time, na.rm = TRUE))

# Produce table with p value comparing median times 
wilcox_median_time <- Pack_year_cessation_bins_combined %>%
  wilcox_test(Time ~ Pack_year_bins) %>%
  adjust_pvalue(method = "BH") %>%
  arrange(p)

# Figure 5 ######
# Figure 5A
#Threshold smoking pack years for smoking-related tumour in ever-smokers
#SBS4
ever_smoker_SBS4_Pack_year <- data.frame(cbind("SBS4-positive", ever_smoker_SBS4$Smoking..pack.years.))
ever_smoker_no_SBS4_Pack_year <- data.frame(cbind("SBS4-negative", ever_smoker_no_SBS4$Smoking..pack.years.))
ever_smoker_SBS_Pack_year <- rbind(ever_smoker_SBS4_Pack_year, ever_smoker_no_SBS4_Pack_year)
colnames(ever_smoker_SBS_Pack_year) <- c("Signature", "Pack_year")
tapply(as.numeric(ever_smoker_SBS_Pack_year$Pack_year), ever_smoker_SBS_Pack_year$Signature, summary)

#DBS2
ever_smoker_DBS2_Pack_year <- data.frame(cbind("DBS2-positive", ever_smoker_DBS2$Smoking..pack.years.))
ever_smoker_no_DBS2_Pack_year <- data.frame(cbind("DBS2-negative", ever_smoker_no_DBS2$Smoking..pack.years.))
ever_smoker_DBS_Pack_year <- rbind(ever_smoker_DBS2_Pack_year, ever_smoker_no_DBS2_Pack_year)
colnames(ever_smoker_DBS_Pack_year) <- c("Signature", "Pack_year")
tapply(as.numeric(ever_smoker_DBS_Pack_year$Pack_year), ever_smoker_DBS_Pack_year$Signature, summary)

#ID3
ever_smoker_ID3_Pack_year <- data.frame(cbind("ID3-positive", ever_smoker_ID3$Smoking..pack.years.))
ever_smoker_no_ID3_Pack_year <- data.frame(cbind("ID3-negative", ever_smoker_no_ID3$Smoking..pack.years.))
ever_smoker_IDS_Pack_year <- rbind(ever_smoker_ID3_Pack_year, ever_smoker_no_ID3_Pack_year)
colnames(ever_smoker_IDS_Pack_year) <- c("Signature", "Pack_year")
tapply(as.numeric(ever_smoker_IDS_Pack_year$Pack_year), ever_smoker_IDS_Pack_year$Signature, summary)

#Combined smoking signatures and and TMB boxplots
ever_smoker_threshold_df <- data.frame(rbind(ever_smoker_SBS4_Pack_year, ever_smoker_no_SBS4_Pack_year,
                                             ever_smoker_DBS2_Pack_year, ever_smoker_no_DBS2_Pack_year,
                                             ever_smoker_ID3_Pack_year, ever_smoker_no_ID3_Pack_year))
colnames(ever_smoker_threshold_df) <- c("Signatures", "Pack_years")
ever_smoker_threshold_df <- na.omit(ever_smoker_threshold_df)
ever_smoker_threshold_df$Signatures <- as.factor(ever_smoker_threshold_df$Signatures)
ever_smoker_threshold_df$Pack_years <- as.numeric(ever_smoker_threshold_df$Pack_years)
ever_smoker_threshold_df$Signatures <- factor(ever_smoker_threshold_df$Signatures, 
                                              levels = c("SBS4-positive", "SBS4-negative", "DBS2-positive", 
                                                         "DBS2-negative", "ID3-positive", "ID3-negative"))

tapply(ever_smoker_threshold_df$Pack_years, ever_smoker_threshold_df$Signatures, summary)

pack_years_stat <- pairwise.wilcox.test(ever_smoker_threshold_df$Pack_years, ever_smoker_threshold_df$Signatures, p.adjust.method = "BH")

my_comparisons <- split(t(combn(levels(ever_smoker_threshold_df$Signatures), 2)), 
                        seq(nrow(t(combn(levels(ever_smoker_threshold_df$Signatures), 2))))) #all pairwise combinations
my_comparisons <- my_comparisons[c(1,10,15)] #pairwise combinations manual order

# calculate medians for smoking pack years in SBS4 positive or negative groups
ever_smoker_threshold_df %>%
  group_by(Signatures) %>%
  summarise(median_value = median(Pack_years, na.rm = TRUE))

#Pack-years in ever-smokers according to smoking signature status
Smoking_signature_threshold_plot <- ggplot(ever_smoker_threshold_df, aes(x=Signatures, y=Pack_years, fill=Signatures)) + 
  geom_boxplot(colour = "black",  outlier.shape = NA) +
  geom_jitter() +
  theme_classic() +
  labs(y="Smoking pack-years") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 20), legend.position = "none") +
  theme(axis.title.x = element_blank()) +
  theme(axis.text.x = element_text(size = 20)) +
  theme(axis.title.y = element_text(face = "bold", size = 20)) +
  theme(axis.text.y = element_text(size = 20)) +
  scale_fill_manual(values = c("firebrick", "grey", "darkred", "grey", "bisque2", "grey")) +
  scale_y_continuous(limits = c(0,185), breaks = seq(0, 185, by = 20)) +
  geom_signif(comparisons = my_comparisons, textsize=8, fontface = "bold", y_position = 170, vjust = -0.5,
              annotations = c(paste("p =", paste(signif(pack_years_stat$p.value[1,1], digits = 2))),
                              paste("p =", paste(signif(pack_years_stat$p.value[3,3], digits = 2))),
                              paste("p =", paste(signif(pack_years_stat$p.value[5,5], digits = 2))))) #annotations means manually enter the p-values

ggsave2(filename = "Figure_5A.pdf", plot = Smoking_signature_threshold_plot, height = 8, width = 15, units = "in", path = "~/file/path/Figures_final/", device = cairo_pdf)

# Figure 5B
#SBS4 counts according to pack years
Pack_year_SBS_counts <- lung_patients_combined_SBS_counts[,c("Smoking.status", "Smoking..pack.years.","SBS4_counts")]
colnames(Pack_year_SBS_counts) <- c("Smoking_status", "Pack_years", "SBS4_counts")
Pack_year_SBS_counts <- Pack_year_SBS_counts[which(!is.na(Pack_year_SBS_counts$Pack_years)),]

Pack_year_SBS_counts_stat <- summary(lm(Pack_year_SBS_counts$Pack_years ~ Pack_year_SBS_counts$SBS4_counts))

SBS4_counts_Pack_years <- ggplot(data = Pack_year_SBS_counts, aes(x = Pack_years, y = SBS4_counts)) +
  geom_point(aes(colour = Smoking_status), size = 3) +
  theme_classic() +
  labs(x = "Smoking pack-years", y="SBS4 mutations") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 20), legend.position = "right",
        legend.title = element_text(face = "bold", size = 15),legend.text = element_text(size = 15),
        axis.title.x = element_text(face = "bold", size = 20), axis.text.x = element_text(size = 15), 
        axis.title.y = element_text(face = "bold", size = 20), axis.text.y = element_text(size = 15)) +
  scale_x_continuous(limits = c(0,165), breaks = seq(0, 165, by=20)) +
  scale_y_continuous(limits = c(0,160000), breaks = seq(0, 160000, by=20000)) +
  scale_colour_manual(values = c("Current" = "firebrick", "Ex" = "steelblue", 
                                 "Never" = "forestgreen"), name = "Smoking status") +
  geom_vline(data = Pack_year_SBS_counts, color = "red", linetype = "longdash", aes(xintercept = median(Pack_years, na.rm = T))) +
  geom_smooth(method = lm, se = T, formula = y ~ x) + #linear regression 
  annotate("text", x = 140, y = 157500, size = 6, label = paste("R\u00b2 adj =", paste(signif(Pack_year_SBS_counts_stat$adj.r.squared, digits = 2)))) +
  annotate("text", x = 140, y = 150000, size = 6, label = paste("p =", paste(signif(Pack_year_SBS_counts_stat$coefficients[2,4], digits = 2))))

ggsave2(filename = "Figure_5B.pdf", plot = SBS4_counts_Pack_years, height = 8, width = 8, units = "in", path = "~/file/path/Figures_final/", device = cairo_pdf)

#Figure 5C
#DBS2 counts according to pack years
Pack_year_DBS_counts <- lung_patients_combined_DBS_counts[,c("Smoking.status", "Smoking..pack.years.","DBS2_counts")]
colnames(Pack_year_DBS_counts) <- c("Smoking_status", "Pack_years", "DBS2_counts")
Pack_year_DBS_counts <- Pack_year_DBS_counts[which(!is.na(Pack_year_DBS_counts$Pack_years)),]

Pack_year_DBS_counts_stat <- summary(lm(Pack_year_DBS_counts$Pack_years ~ Pack_year_DBS_counts$DBS2_counts))

DBS2_counts_Pack_years <- ggplot(data = Pack_year_DBS_counts, aes(x = Pack_years, y = DBS2_counts)) +
  geom_point(aes(colour = Smoking_status), size = 3) +
  theme_classic() +
  labs(x = "Smoking pack-years", y="DBS2 mutations") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 20), legend.position = "right",
        legend.title = element_text(face = "bold", size = 15),legend.text = element_text(size = 15),
        axis.title.x = element_text(face = "bold", size = 20), axis.text.x = element_text(size = 15), 
        axis.title.y = element_text(face = "bold", size = 20), axis.text.y = element_text(size = 15)) +
  scale_x_continuous(limits = c(0,170), breaks = seq(0, 170, by=20)) +
  scale_y_continuous(limits = c(0,4500), breaks = seq(0, 4500, by=500)) +
  scale_colour_manual(values = c("Current" = "firebrick", "Ex" = "steelblue", 
                                 "Never" = "forestgreen"), name = "Smoking status") +
  geom_vline(data = Pack_year_DBS_counts, color = "red", linetype = "longdash", aes(xintercept = median(Pack_years, na.rm = T))) +
  geom_smooth(method = lm, se = T, formula = y ~ x) + #linear regression 
  annotate("text", x = 140, y = 4500, size = 6, label = paste("R\u00b2 adj =", paste(signif(Pack_year_DBS_counts_stat$adj.r.squared, digits = 2)))) +
  annotate("text", x = 140, y = 4300, size = 6, label = paste("p =", paste(signif(Pack_year_DBS_counts_stat$coefficients[2,4], digits = 2))))

ggsave2(filename = "Figure_5C.pdf", plot = DBS2_counts_Pack_years, height = 8, width = 8, units = "in", path = "~/file/path/Figures_final/", device = cairo_pdf)

#Figure 5D
#ID3 counts according to pack years
Pack_year_IDS_counts <- lung_patients_combined_IDS_counts[,c("Smoking.status", "Smoking..pack.years.","ID3_counts")]
colnames(Pack_year_IDS_counts) <- c("Smoking_status", "Pack_years", "ID3_counts")
Pack_year_IDS_counts <- Pack_year_IDS_counts[which(!is.na(Pack_year_IDS_counts$Pack_years)),]

Pack_year_IDS_counts_stat <- summary(lm(Pack_year_IDS_counts$Pack_years ~ Pack_year_IDS_counts$ID3_counts))

ID3_counts_Pack_years <- ggplot(data = Pack_year_IDS_counts, aes(x = Pack_years, y = ID3_counts)) +
  geom_point(aes(colour = Smoking_status), size = 3) +
  theme_classic() +
  labs(x = "Smoking pack-years", y="ID3 mutations") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 20), legend.position = "right",
        legend.title = element_text(face = "bold", size = 15),legend.text = element_text(size = 15),
        axis.title.x = element_text(face = "bold", size = 20), axis.text.x = element_text(size = 15), 
        axis.title.y = element_text(face = "bold", size = 20), axis.text.y = element_text(size = 15)) +
  scale_x_continuous(limits = c(0,170), breaks = seq(0, 170, by=20)) +
  scale_y_continuous(limits = c(0,8000), breaks = seq(0, 8000, by=1000)) +
  scale_colour_manual(values = c("Current" = "firebrick", "Ex" = "steelblue", 
                                 "Never" = "forestgreen"), name = "Smoking status") +
  geom_vline(data = Pack_year_IDS_counts, color = "red", linetype = "longdash", aes(xintercept = median(Pack_years, na.rm = T))) +
  geom_smooth(method = lm, se = T, formula = y ~ x) + #linear regression 
  annotate("text", x = 140, y = 7800, size = 6, label = paste("R\u00b2 adj =", paste(signif(Pack_year_IDS_counts_stat$adj.r.squared, digits = 2)))) +
  annotate("text", x = 140, y = 7500, size = 6, label = paste("p =", paste(signif(Pack_year_IDS_counts_stat$coefficients[2,4], digits = 2))))

ggsave2(filename = "Figure_5D.pdf", plot = ID3_counts_Pack_years, height = 8, width = 8, units = "in", path =  "~/file/path/Figures_final/", device = cairo_pdf)

combined_pack_years_plot <- plot_grid(SBS4_counts_Pack_years,DBS2_counts_Pack_years,ID3_counts_Pack_years, align = "h", nrow = 1)

ggsave2(filename = "Figure_5BCD_pdf", plot = combined_pack_years_plot, height = 8, width = 22, units = "in", path = "~/file/path/Figures_final/", device = cairo_pdf)

# Supp Figure 4 #####
# SBS4 activity corrected for age 
# Figure S4A
Pack_year_age_SBS4 <- lung_patients_combined_SBS_counts[,c("Smoking.status","Smoking..pack.years.","SBS4_counts", "Age.at.diagnosis")]
colnames(Pack_year_age_SBS4) <- c("Smoking_status", "Pack_year", "SBS4", "Age")
Pack_year_age_SBS4 <- Pack_year_age_SBS4[!is.na(Pack_year_age_SBS4$Pack_year),]
Pack_year_age_SBS4$Pack_year_corrected <- Pack_year_age_SBS4$Pack_year/Pack_year_age_SBS4$Age

summary(Pack_year_age_SBS4$Pack_year_corrected)
Pack_year_SBS4_stat <- summary(lm(Pack_year_age_SBS4$Pack_year_corrected ~ Pack_year_age_SBS4$SBS4))

SBS4_counts_Pack_years <- ggplot(data = Pack_year_age_SBS4, aes(x = Pack_year_corrected, y = SBS4)) +
  geom_point(aes(colour = Smoking_status), size = 3) +
  theme_classic() +
  labs(x = "Pack-years corrected for age", y="SBS4 mutations") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 20), legend.position = "right",
        legend.title = element_text(face = "bold", size = 15),legend.text = element_text(size = 15),
        axis.title.x = element_text(face = "bold", size = 20), axis.text.x = element_text(size = 15), 
        axis.title.y = element_text(face = "bold", size = 20), axis.text.y = element_text(size = 15)) +
  scale_x_continuous(limits = c(0,2.4), breaks = seq(0, 2.4, by=0.3)) +
  scale_y_continuous(limits = c(0,160000), breaks = seq(0, 160000, by=20000)) +
  scale_colour_manual(values = c("Current" = "firebrick", "Ex" = "steelblue", 
                                 "Never" = "forestgreen"), name = "Smoking status") +
  geom_vline(data = Pack_year_age_SBS4, color = "red", linetype = "longdash", aes(xintercept = median(Pack_year_age_SBS4$Pack_year_corrected, na.rm = T))) +
  geom_smooth(method = lm, se = T, formula = y ~ x) + #linear regression 
  annotate("text", x = 2, y = 157500, size = 6, label = paste("R\u00b2 adj =", paste(signif(Pack_year_SBS4_stat$adj.r.squared, digits = 2)))) +
  annotate("text", x = 2, y = 150000, size = 6, label = paste("p =", paste(signif(Pack_year_SBS4_stat$coefficients[2,4], digits = 2))))

ggsave2(filename = "Figure_S4A.pdf", plot = SBS4_counts_Pack_years, height = 8, width = 8, units = "in", path = "~/file/path/Figures_final/", device = cairo_pdf)


# Figure S4B
Pack_year_age_DBS2 <- lung_patients_combined_DBS_counts[,c("Smoking.status","Smoking..pack.years.","DBS2_counts", "Age.at.diagnosis")]
colnames(Pack_year_age_DBS2) <- c("Smoking_status", "Pack_year", "DBS2", "Age")
Pack_year_age_DBS2 <- Pack_year_age_DBS2[!is.na(Pack_year_age_DBS2$Pack_year),]
Pack_year_age_DBS2$Pack_year_corrected <- Pack_year_age_DBS2$Pack_year/Pack_year_age_DBS2$Age

summary(Pack_year_age_DBS2$Pack_year_corrected)
Pack_year_DBS2_stat <- summary(lm(Pack_year_age_DBS2$Pack_year_corrected ~ Pack_year_age_DBS2$DBS2))

DBS2_counts_Pack_years <- ggplot(data = Pack_year_age_DBS2, aes(x = Pack_year_corrected, y = DBS2)) +
  geom_point(aes(colour = Smoking_status), size = 3) +
  theme_classic() +
  labs(x = "Pack-years corrected for age", y="DBS2 mutations") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 20), legend.position = "right",
        legend.title = element_text(face = "bold", size = 15),legend.text = element_text(size = 15),
        axis.title.x = element_text(face = "bold", size = 20), axis.text.x = element_text(size = 15), 
        axis.title.y = element_text(face = "bold", size = 20), axis.text.y = element_text(size = 15)) +
  scale_x_continuous(limits = c(0,2.4), breaks = seq(0, 2.4, by=0.3)) +
  scale_y_continuous(limits = c(0,4500), breaks = seq(0, 4500, by=500)) +
  scale_colour_manual(values = c("Current" = "firebrick", "Ex" = "steelblue", 
                                 "Never" = "forestgreen"), name = "Smoking status") +
  geom_vline(data = Pack_year_age_DBS2, color = "red", linetype = "longdash", aes(xintercept = median(Pack_year_age_DBS2$Pack_year_corrected, na.rm = T))) +
  geom_smooth(method = lm, se = T, formula = y ~ x) + #linear regression 
  annotate("text", x = 2, y = 4500, size = 6, label = paste("R\u00b2 adj =", paste(signif(Pack_year_DBS2_stat$adj.r.squared, digits = 2)))) +
  annotate("text", x = 2, y = 4300, size = 6, label = paste("p =", paste(signif(Pack_year_DBS2_stat$coefficients[2,4], digits = 2))))

ggsave2(filename = "Figure_S4B.pdf", plot = DBS2_counts_Pack_years, height = 8, width = 8, units = "in", path = "~/file/path/Figures_final/", device = cairo_pdf)

# Figure S4C
Pack_year_age_ID3 <- lung_patients_combined_IDS_counts[,c("Smoking.status","Smoking..pack.years.","ID3_counts", "Age.at.diagnosis")]
colnames(Pack_year_age_ID3) <- c("Smoking_status", "Pack_year", "ID3", "Age")
Pack_year_age_ID3 <- Pack_year_age_ID3[!is.na(Pack_year_age_ID3$Pack_year),]
Pack_year_age_ID3$Pack_year_corrected <- Pack_year_age_ID3$Pack_year/Pack_year_age_ID3$Age

summary(Pack_year_age_ID3$Pack_year_corrected)
Pack_year_ID3_stat <- summary(lm(Pack_year_age_ID3$Pack_year_corrected ~ Pack_year_age_ID3$ID3))

ID3_counts_Pack_years <- ggplot(data = Pack_year_age_ID3, aes(x = Pack_year_corrected, y = ID3)) +
  geom_point(aes(colour = Smoking_status), size = 3) +
  theme_classic() +
  labs(x = "Pack-years corrected for age", y="ID3 mutations") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 20), legend.position = "right",
        legend.title = element_text(face = "bold", size = 15),legend.text = element_text(size = 15),
        axis.title.x = element_text(face = "bold", size = 20), axis.text.x = element_text(size = 15), 
        axis.title.y = element_text(face = "bold", size = 20), axis.text.y = element_text(size = 15)) +
  scale_x_continuous(limits = c(0,2.4), breaks = seq(0, 2.4, by=0.3)) +
  scale_y_continuous(limits = c(0,8000), breaks = seq(0, 8000, by=1000)) +
  scale_colour_manual(values = c("Current" = "firebrick", "Ex" = "steelblue", 
                                 "Never" = "forestgreen"), name = "Smoking status") +
  geom_vline(data = Pack_year_age_ID3, color = "red", linetype = "longdash", aes(xintercept = median(Pack_year_age_ID3$Pack_year_corrected, na.rm = T))) +
  geom_smooth(method = lm, se = T, formula = y ~ x) + #linear regression 
  annotate("text", x = 2, y = 7800, size = 6, label = paste("R\u00b2 adj =", paste(signif(Pack_year_ID3_stat$adj.r.squared, digits = 2)))) +
  annotate("text", x = 2, y = 7500, size = 6, label = paste("p =", paste(signif(Pack_year_ID3_stat$coefficients[2,4], digits = 2))))

ggsave2(filename = "Figure_S4C.pdf", plot = ID3_counts_Pack_years, height = 8, width = 8, units = "in", path = "~/file/path/Figures_final/", device = cairo_pdf)

combined_pack_years_corrected_plot <- plot_grid(SBS4_counts_Pack_years,DBS2_counts_Pack_years,ID3_counts_Pack_years, align = "h", nrow = 1)

ggsave2(filename = "Figure_4FtoG_corrected.pdf", plot = combined_pack_years_corrected_plot, height = 8, width = 22, units = "in", path = "~/file/path/Figures_final/", device = cairo_pdf)

# Supp 4D 
Pack_year_TMB <- Pack_year_age[,c("Smoking.status", "Smoking..pack.years.","TMB")]
colnames(Pack_year_TMB) <- c("Smoking_status", "Pack_year", "TMB")
Pack_year_TMB <- Pack_year_TMB[!is.na(Pack_year_TMB$Pack_year),]

Pack_year_TMB_stat <- summary(lm(Pack_year_TMB$Pack_year ~ Pack_year_TMB$TMB))

TMB_Pack_years_plot <- ggplot(data = Pack_year_TMB, aes(x = Pack_year, y = TMB)) +
  geom_point(aes(colour = Smoking_status), size = 3) +
  theme_classic() +
  labs(x = "Smoking pack-years", y="TMB (mutations/Mb)") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 20), legend.position = "right",
        legend.title = element_text(face = "bold", size = 15),legend.text = element_text(size = 15),
        axis.title.x = element_text(face = "bold", size = 20), axis.text.x = element_text(size = 15), 
        axis.title.y = element_text(face = "bold", size = 20), axis.text.y = element_text(size = 15)) +
  scale_x_continuous(limits = c(0,170), breaks = seq(0, 170, by=20)) +
  scale_y_continuous(limits = c(0,60), breaks = seq(0, 60, by=5)) +
  scale_colour_manual(values = c("Current" = "firebrick", "Ex" = "steelblue", 
                                 "Never" = "forestgreen"), name = "Smoking status") +
  geom_vline(data = Pack_year_TMB, color = "red", linetype = "longdash", aes(xintercept = median(Pack_year_TMB$Pack_year, na.rm = T))) +
  geom_smooth(method = lm, se = T, formula = 'y ~ x') + #linear regression 
  annotate("text", x = 150, y = 50, size = 6, label = paste("R\u00b2 adj =", paste(signif(Pack_year_TMB_stat$adj.r.squared, digits = 2)))) +
  annotate("text", x = 150, y = 48, size = 6, label = paste("p =", paste(signif(Pack_year_TMB_stat$coefficients[2,4], digits = 2))))

ggsave2(filename = "Figure_S4D.pdf", plot = TMB_Pack_years_plot, height = 8, width = 8, units = "in", path = "~/file/path/Figures_final/", device = cairo_pdf)

# Figure Supp 4E pack years corrected for by age
Pack_year_age <- lung_patients_combined
Pack_year_age_TMB <- Pack_year_age[,c("Smoking..pack.years.","TMB", "Age.at.diagnosis")]
colnames(Pack_year_age_TMB) <- c("Pack_year", "TMB", "Age")
Pack_year_age_TMB <- Pack_year_age_TMB[!is.na(Pack_year_age_TMB$Pack_year),]
Pack_year_age_TMB$Pack_year_corrected <- Pack_year_age_TMB$Pack_year/Pack_year_age_TMB$Age

summary(Pack_year_age_TMB$Pack_year_corrected)
Pack_year_TMB_stat_corrected <- summary(lm(Pack_year_age_TMB$Pack_year_corrected ~ Pack_year_age_TMB$TMB))

Pack_year_age_TMB_corrected_plot <- ggplot(data = Pack_year_age_TMB, aes(x = Pack_year_corrected, y = TMB)) +
  geom_point() +
  theme_classic() +
  labs(x = "Pack-years corrected for age", y="TMB (mutations/Mb)") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 20), legend.position = "none") +
  theme(axis.title.x = element_text(face = "bold", size = 20)) +
  theme(axis.text.x = element_text(size = 20)) +
  theme(axis.title.y = element_text(face = "bold", size = 20)) +
  theme(axis.text.y = element_text(size = 20)) +
  scale_x_continuous(limits = c(0,2.4), breaks = seq(0, 2.4, by=0.3)) +
  scale_y_continuous(limits = c(0,50), breaks = seq(0, 50, by=5)) +
  geom_vline(data = Pack_year_age_TMB, color = "red", linetype = "longdash", aes(xintercept = median(Pack_year_age_TMB$Pack_year_corrected, na.rm = T))) +
  geom_smooth(method = lm, se = T, formula = y ~ x) + #linear regression 
  annotate("text", x = 2, y = 50, size = 8, label = paste("R\u00b2 adj =", paste(signif(Pack_year_TMB_stat_corrected$adj.r.squared, digits = 2)))) +
  annotate("text", x = 2, y = 48, size = 8, label = paste("p =", paste(signif(Pack_year_TMB_stat_corrected$coefficients[2,4], digits = 2))))

ggsave2(filename = "Figure_S4E.pdf", plot = Pack_year_age_TMB_corrected_plot, height = 8, width = 8, units = "in", path = "~/file/path/Figures_final/", device = cairo_pdf)
