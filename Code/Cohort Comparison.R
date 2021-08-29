## Microarray cohort and single cell cohort comparison

library(tidyverse)
library(ggplot2)
library(rstatix)
library(ggpubr)
library(scales)

setwd("~/Cohort Comparison/")

#Prepare dataset
micro_pheno <- read.csv("~/VascMicroarray/Usable Data/Phenodata.csv")
micro_data <- micro_pheno %>% 
  subset(cohort == "Vasculitis") %>% 
  dplyr::select(study_id, sex, age, timepoint_short, BVAS, crp, gfr, 
         pr3_titre, ent, renal, joint)
cohort <- replicate(nrow(micro_data), "microarray")
micro_data <- cbind(cohort, micro_data)


single_pheno <- read.csv("~/Single_Cell/single_cell_phenodata.csv")
single_data <- single_pheno %>% 
  subset(tissue == "PBMC") %>% 
  subset(timepoint != "control") %>% 
  dplyr::select(Trial_ID, sex, age, timepoint, bvas_score, crp, gfr, 
                pr3_titre, ent, renal, joints)
#Change timepoint names to match micro data
single_data$timepoint <- str_replace_all(single_data$timepoint, 'D1', 'zero')
single_data$timepoint <- str_replace_all(single_data$timepoint, 'W12', 'three')
cohort <- replicate(nrow(single_data), "singlecell")
single_data <- cbind(cohort, single_data)
colnames(single_data) <- colnames(micro_data)

all_data <- rbind(micro_data, single_data)
str(all_data)
all_data <- all_data %>% 
  transform(age = as.numeric(age),
            BVAS = as.numeric(BVAS),
            crp = as.numeric(crp),
            gfr = as.numeric(gfr),
            pr3_titre = as.numeric(pr3_titre),
            ent = as.numeric(ent),
            renal = as.numeric(renal),
            joint = as.numeric(joint)) 
all_data$cohort <- factor(all_data$cohort, levels = c("microarray", "singlecell"))

  

#All AAV patients - t_test
t_age <- all_data %>% 
  t_test(age ~ cohort)

t_bvas <- all_data %>% 
  t_test(BVAS ~ cohort)

t_crp <- all_data %>% 
  t_test(crp ~ cohort)

t_gfr <- all_data %>% 
  t_test(gfr ~ cohort)

t_pr3_titre <- all_data %>% 
  t_test(pr3_titre ~ cohort)

all_t <- rbind(t_age, t_bvas, t_crp, t_gfr, t_pr3_titre)

#Active vs active and remission vs remission
time_data <- subset(all_data, timepoint_short != "twelve")

time_t_age <- time_data %>% 
  group_by(timepoint_short) %>% 
  t_test(age ~ cohort) %>% 
  ungroup()

time_t_bvas <- time_data %>% 
  group_by(timepoint_short) %>% 
  t_test(BVAS ~ cohort) %>% 
  ungroup()

time_t_crp <- time_data %>% 
  group_by(timepoint_short) %>% 
  t_test(crp ~ cohort) %>% 
  ungroup()

time_t_gfr <- time_data %>% 
  group_by(timepoint_short) %>% 
  t_test(gfr ~ cohort) %>% 
  ungroup()

time_t_pr3_titre <- time_data %>% 
  group_by(timepoint_short) %>% 
  t_test(pr3_titre ~ cohort) %>% 
  ungroup()

time_t <- rbind(time_t_age, time_t_bvas, time_t_crp, time_t_gfr, time_t_pr3_titre)


##Comparing categorical data in all AAV
#Create dataframe for sex
sex <- table(all_data$cohort, all_data$sex)
dimnames(sex)
f_sex <- fisher_test(sex, detailed = TRUE)

ent <- table(all_data$cohort, all_data$ent)
dimnames(ent)
f_ent <- fisher_test(ent, detailed = TRUE)

renal <- table(all_data$cohort, all_data$renal)
dimnames(renal)
f_renal <- fisher_test(renal, detailed = TRUE)

joint <- table(all_data$cohort, all_data$joint)
dimnames(joint)
f_joint <- fisher_test(joint, detailed = TRUE)


all_f <- rbind(f_sex, f_ent, f_renal, f_joint)
var <- c("sex", "ent", "renal", "joint")
all_f <- cbind(var, all_f)

#In Active AAV
active_data <- subset(all_data, timepoint_short == "zero")

active_sex <- table(active_data$cohort, active_data$sex)
dimnames(active_sex)
active_f_sex <- fisher_test(active_sex, detailed = TRUE)

active_ent <- table(active_data$cohort, active_data$ent)
dimnames(active_ent)
active_f_ent <- fisher_test(active_ent, detailed = TRUE)

active_renal <- table(active_data$cohort, active_data$renal)
dimnames(active_renal)
active_f_renal <- fisher_test(active_renal, detailed = TRUE)

active_joint <- table(active_data$cohort, active_data$joint)
dimnames(active_joint)
active_f_joint <- fisher_test(active_joint, detailed = TRUE)


active_f <- rbind(active_f_sex, active_f_ent, active_f_renal, active_f_joint)
active_f <- cbind(var, active_f)


#In Remis AAV
remis_data <- subset(all_data, timepoint_short == "three")

remis_sex <- table(remis_data$cohort, remis_data$sex)
dimnames(remis_sex)
remis_f_sex <- fisher_test(remis_sex, detailed = TRUE)

remis_ent <- table(remis_data$cohort, remis_data$ent)
dimnames(remis_ent)
remis_f_ent <- fisher_test(remis_ent, detailed = TRUE)

remis_renal <- table(remis_data$cohort, remis_data$renal)
dimnames(remis_renal)
remis_f_renal <- fisher_test(remis_renal, detailed = TRUE)

remis_joint <- table(remis_data$cohort, remis_data$joint)
dimnames(remis_joint)
remis_f_joint <- fisher_test(remis_joint, detailed = TRUE)


remis_f <- rbind(remis_f_sex, remis_f_ent, remis_f_renal, remis_f_joint)
remis_f <- cbind(var, remis_f)


  