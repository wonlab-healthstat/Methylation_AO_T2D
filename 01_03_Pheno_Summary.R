#####################################
### Study : Bohun-Methylation
### Purpose of program : Baseline demographics summary
### Pre-required program : NA
### Output files : Baseline_KARE.csv
### Programmed by : Sujin Seo
### Draft Date : 2022-06-03
### Revision : NA
### Software (version) : R (v3.6.3)
#####################################

#####################################
##### Load Required Libraries
#####################################
library(dplyr)
library(readxl)

#####################################
##### Set Working Directories
#####################################

dr_dir <- "/home2/sjseo/Bohun/Methylation/derived_data/"

#####################################
##### Read data
#####################################

anal_dat <- read.table(paste0(dr_dir, "03_02_QC_Pheno_RE.csv"), sep = ",", header = T) %>%
  mutate(DM_ANAL = DM + 1)
dat <- read.table(paste0(dr_dir, "01_02_Baseline_KARE.csv"), sep = ",", header = T) %>%
  left_join(., anal_dat %>% select(Sample_ID, DM_ANAL), by = c("DIST_ID" = "Sample_ID")) %>%
  mutate(DM_IDEN = ifelse(DM == DM_ANAL, "T", "F"),
         cohort = "KARE", ID = DIST_ID, Sex = AS5_SEX,
         DM = DM_ANAL)
 
bohun <- read_xlsx("/data/bohun_methylation/EPC array total information_rex.xlsx") %>%
  mutate(cohort = "BH", ID = Sample_ID, DM = DM+1, Sex = 1)

#####################################
##### Summarize
#####################################
[1] "DIST_ID"        "AS5_AGE"        "AS5_SEX"        "AS5_BMI"        "AS5_EDATE5"     "AS5_HBA1C"      "AS5_CREATININE" "AS5_eGRF"       "BIRTH"          "DM"             "DMDIAGDATE"    
[12] "DMONSET"        "DMDUR"          "MI"             "MICU"           "HTN"            "HTNCU"          "LCA"            "LCACU"          "GCA"            "GCACU"          "HCC"           
[23] "HCCCU"          "COLCA"          "COLCACU"        "PACA"           "PACACU"         "CEVA"           "CEVACU"         "DM_ANAL"        "DM_IDEN"    


> colnames(bohun)
[1] "Sample_ID"   "id"          "name"        "sex"         "age"         "group"       "DM"          "DM_duration" "DM_age"      "HbA1c"       "eGFR"        "cr"          "BMI"        
[14] "AMI"         "ca"          "HTN"         "PD"          "Stoke"       "Liver"       "skin"        "renal"      


##### Age
a <- bohun %>% select(cohort, ID, DM, Sex, Age = age)
b <- dat %>% select(cohort, ID, DM, Sex, Age = AS5_AGE)
c <- rbind(a, b) %>% mutate(Group = ifelse(cohort == "BH", "A", 
                                           ifelse(DM == 2, "B", "C")))

c %>%  summarize(n = length(unique(ID)), m = format(mean(Age), nsmall = 2), sd = format(sd(Age), nsmall = 2))
c %>% group_by(cohort, DM, Sex) %>% summarize(n = length(unique(ID)), m = format(mean(Age), nsmall = 2), sd = format(sd(Age), nsmall = 2))
summary(aov(Age ~ Group, data = c))
TukeyHSD(aov(Age ~ Group, data = c))

##### BMI
a <- bohun %>% select(cohort, ID, DM, Sex, BMI)
b <- dat %>% select(cohort, ID, DM, Sex, BMI = AS5_BMI)
c <- rbind(a, b) %>% mutate(Group = ifelse(cohort == "BH", "A", 
                                           ifelse(DM == 2, "B", "C")))

c %>%  summarize(n = length(unique(ID)), m = format(mean(BMI), nsmall = 2), sd = format(sd(BMI), nsmall = 2))
c %>% group_by(cohort, DM, Sex) %>% summarize(n = length(unique(ID)), m = format(mean(BMI), nsmall = 2), sd = format(sd(BMI), nsmall = 2))
summary(aov(BMI ~ Group, data = c))
TukeyHSD(aov(BMI ~ Group, data = c))


##### Duration DM
a <- bohun %>% select(cohort, ID, DM, Sex, DMDUR = DM_duration)
b <- dat %>% select(cohort, ID, DM, Sex, DMDUR)
c <- rbind(a, b) %>% mutate(Group = ifelse(cohort == "BH", "A", 
                                           ifelse(DM == 2, "B", "C")))

c %>%  summarize(n = length(unique(ID)), m = format(mean(DMDUR, na.rm = T), nsmall = 2), sd = format(sd(DMDUR, na.rm = T), nsmall = 2))
c %>% group_by(cohort, DM, Sex) %>% 
  summarize(n = length(unique(ID)), m = format(mean(DMDUR, na.rm = T), nsmall = 2), sd = format(sd(DMDUR, na.rm = T), nsmall = 2))

summary(aov(DMDUR ~ Group, data = c))
t.test(a$DMDUR, b$DMDUR)

##### Onset of DM
a <- bohun %>% select(cohort, ID, DM, Sex, DMONSET = DM_age)
b <- dat %>% select(cohort, ID, DM, Sex, DMONSET)
c <- rbind(a, b) %>% mutate(Group = ifelse(cohort == "BH", "A", 
                                           ifelse(DM == 2, "B", "C")))

c %>%  summarize(n = length(unique(ID)), m = format(mean(DMONSET, na.rm = T), nsmall = 2), sd = format(sd(DMONSET, na.rm = T), nsmall = 2))
c %>% group_by(cohort, DM, Sex) %>% 
  summarize(n = length(unique(ID)), m = format(mean(DMONSET, na.rm = T), nsmall = 2), sd = format(sd(DMONSET, na.rm = T), nsmall = 2))
summary(aov(DMONSET ~ Group, data = c))
t.test(a$DMONSET, b$DMONSET)

##### HbA1c
a <- bohun %>% select(cohort, ID, DM, Sex, HBA1C = HbA1c)
b <- dat %>% select(cohort, ID, DM, Sex, HBA1C = AS5_HBA1C)
c <- rbind(a, b) %>% mutate(Group = ifelse(cohort == "BH", "A", 
                                           ifelse(DM == 2, "B", "C")))

c %>%  summarize(n = length(unique(ID)), m = format(mean(HBA1C, na.rm = T), nsmall = 2), sd = format(sd(HBA1C, na.rm = T), nsmall = 2))
c %>% group_by(cohort, DM, Sex) %>% 
  summarize(n = length(unique(ID)), m = format(mean(HBA1C, na.rm = T), nsmall = 2), sd = format(sd(HBA1C, na.rm = T), nsmall = 2))
summary(aov(HBA1C ~ Group, data = c))
plot(TukeyHSD(aov(HBA1C ~ Group, data = c)))

##### eGFR
a <- bohun %>% select(cohort, ID, DM, Sex, eGFR)
b <- dat %>% select(cohort, ID, DM, Sex, eGFR = AS5_eGRF)
c <- rbind(a, b) %>% mutate(Group = ifelse(cohort == "BH", "A", 
                                           ifelse(DM == 2, "B", "C")))

c %>%  summarize(n = length(unique(ID)), m = format(mean(eGFR, na.rm = T), nsmall = 2), sd = format(sd(eGFR, na.rm = T), nsmall = 2))
c %>% group_by(cohort, DM, Sex) %>% 
  summarize(n = length(unique(ID)), m = format(mean(eGFR, na.rm = T), nsmall = 2), sd = format(sd(eGFR, na.rm = T), nsmall = 2))
summary(aov(eGFR ~ Group, data = c))
TukeyHSD(aov(eGFR ~ Group, data = c))

##### Creatinine
a <- bohun %>% select(cohort, ID, DM, Sex, CREAT = cr)
b <- dat %>% select(cohort, ID, DM, Sex, CREAT = AS5_CREATININE)
c <- rbind(a, b) %>% mutate(Group = ifelse(cohort == "BH", "A", 
                                           ifelse(DM == 2, "B", "C")))

c %>%  summarize(n = length(unique(ID)), m = format(mean(CREAT, na.rm = T), nsmall = 2), sd = format(sd(CREAT, na.rm = T), nsmall = 2))
c %>% group_by(cohort, DM, Sex) %>% 
  summarize(n = length(unique(ID)), m = format(mean(CREAT, na.rm = T), nsmall = 2), sd = format(sd(CREAT, na.rm = T), nsmall = 2))
summary(aov(CREAT ~ Group, data = c))
TukeyHSD(aov(CREAT ~ Group, data = c))

##### Comobidity
### MI
a <- bohun %>% select(cohort, ID, DM, Sex, AMI)
b <- dat %>% mutate(AMI = MI-1) %>% select(cohort, ID, DM, Sex, AMI)
c <- rbind(a, b) %>% mutate(Group = ifelse(cohort == "BH", "A", 
                                           ifelse(DM == 2, "B", "C")))

c %>%  summarize(N = length(unique(ID)), n = sum(AMI), r = format(sum(AMI)/length(unique(ID))*100, nsmall = 2))
c %>% group_by(cohort, DM, Sex) %>% 
  summarize(N = length(unique(ID)), n = sum(AMI), r = format(sum(AMI)/length(unique(ID))*100, nsmall = 2))
with(c, table(Group, AMI))
chisq.test(with(c, table(Group, AMI)))
fisher.test(with(c, table(Group, AMI)))

### Hypertension
a <- bohun %>% select(cohort, ID, DM, Sex, HTN)
b <- dat %>% mutate(HTN = HTN-1) %>% select(cohort, ID, DM, Sex, HTN)
c <- rbind(a, b) %>% mutate(Group = ifelse(cohort == "BH", "A", 
                                           ifelse(DM == 2, "B", "C")))

c %>%  summarize(N = length(unique(ID)), n = sum(HTN), r = format(sum(HTN)/length(unique(ID))*100, nsmall = 2))
c %>% group_by(cohort, DM, Sex) %>% 
  summarize(N = length(unique(ID)), n = sum(HTN), r = format(sum(HTN)/length(unique(ID))*100, nsmall = 2))
with(c, table(Group, HTN))
chisq.test(with(c, table(Group, HTN)))
fisher.test(with(c, table(Group, HTN)))

### Cancer
a <- bohun %>% select(cohort, ID, DM, Sex, CANCER = ca)
b <- dat %>% mutate(CANCER = ifelse(sum(LCA-1, GCA-1, HCC-1, COLCA-1, PACA-1)>0, 1, 0)) %>% 
  select(cohort, ID, DM, Sex, CANCER)
c <- rbind(a, b) %>% mutate(Group = ifelse(cohort == "BH", "A", 
                                           ifelse(DM == 2, "B", "C")))

with(c, table(Group, CANCER))
c %>%  summarize(N = length(unique(ID)), n = sum(CANCER), r = format(sum(CANCER)/length(unique(ID))*100, nsmall = 2))
c %>% group_by(cohort, DM, Sex) %>% 
  summarize(N = length(unique(ID)), n = sum(CANCER), r = format(sum(CANCER)/length(unique(ID))*100, nsmall = 2))
chisq.test(with(c, table(Group, CANCER)))
fisher.test(with(c, table(Group, CANCER)))

### Stroke
a <- bohun %>% select(cohort, ID, DM, Sex, CEVA = Stoke)
b <- dat %>% mutate(CEVA = CEVA-1) %>% select(cohort, ID, DM, Sex, CEVA)
c <- rbind(a, b) %>% mutate(Group = ifelse(cohort == "BH", "A", 
                                           ifelse(DM == 2, "B", "C")))

with(c, table(Group, CEVA))
c %>%  summarize(N = length(unique(ID)), n = sum(CEVA), r = format(sum(CEVA)/length(unique(ID))*100, nsmall = 2))
c %>% group_by(cohort, DM, Sex) %>% 
  summarize(N = length(unique(ID)), n = sum(CEVA), r = format(sum(CEVA)/length(unique(ID))*100, nsmall = 2))
chisq.test(with(c, table(Group, CEVA)))
fisher.test(with(c, table(Group, CEVA)))
