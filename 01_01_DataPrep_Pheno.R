#####################################
### Study : Bohun-Methylation
### Purpose of program : Summarize the pheno data
###                     (Cohort, DM, SEX)
### Pre-required program : NA
### Output files : DataPrep_Pheno.csv
### Programmed by : Sujin Seo
### Draft Date : 2021-09-29
### Revision : NA
### Software (version) : R (v3.6.3)
#####################################

#####################################
##### Load Required Libraries
#####################################
library(readxl)
library(writexl)
library(dplyr)
library(readr)
library(plyr)
library(purrr)


#####################################
##### Set Working Directories
#####################################

##### KCDC Codebook
cb_path <- "/data/bohun_methylation/5566/역학정보/3.코드북/"

##### KCDC Pheno
pheno_path <- "/data/bohun_methylation/5566/역학정보/2.EXCEL/지역사회/"
pheno_path2 <- "/data/bohun_methylation/5566/메틸레이션질본_정상_추가"

#### Bohun Pheno
bohun_path <- "/data/bohun_methylation/"

##### Output Directories
out_path <- "/home2/sjseo/Bohun/Methylation/derived_data/"


#############################################
##### KCDC
############################################
dm_var <- c("AS1_PDDM", "AS2_PDDM", "AS3_PDFDM", 
            "AS4_DMDIAG", "AS5_DMDIAG") # Diabetes variable name for each cohort

##### First Data
#pheno_kcdc_1_path <- "C:/Users/SUJIN/Desktop/WonLab/Methylation/Data/KCDC/methlaytion_dat.csv"
pheno_kcdc_1_path <- "/home2/sjseo/Bohun/Methylation/data/KCDC/methlaytion_dat.csv"
pheno_kcdc_1 <- read_csv(pheno_kcdc_1_path, show_col_types = FALSE)

colnames(pheno_kcdc_1) <- toupper(colnames(pheno_kcdc_1))
sex_var <- grep("*(SEX)$", colnames(pheno_kcdc_1), value = TRUE)

pheno_kcdc_1 <- pheno_kcdc_1 %>% 
  select(DIST_ID = MSKID, all_of(c(dm_var, sex_var)))

##### Second Data
### There's no DM information in pheno data so referred to the code books
file_lst <- list.files(cb_path)

extract_dm <- lapply(1:5, function(i) {
  file_nm <- grep(paste0("^\\지역사회기반코호트_\\(", i), file_lst, value = T)
  codebook <- read_excel(paste0(cb_path, file_nm),
                         sheet = "공개코드북", skip = 2)
  
  # DM Information
  tmp_dm_var <- grep(paste0("^\\AS", i), dm_var, value = T)
  tmp_dm_dat <- codebook %>% filter(toupper(`변수명`) == tmp_dm_var) %>% 
    select(DAT = `테이블명\r\n(영문)`, VAR = `변수명`) 
  tmp_dm_file <- paste0(pheno_path, i, "기/", tmp_dm_dat$DAT, ".csv")
  tmp_dm <- read_csv(tmp_dm_file, show_col_types = FALSE) %>% 
    select(DIST_ID, toupper(tmp_dm_dat$VAR))
  
  # Sex Information
  tmp_sex_var <- grep(paste0("^\\AS", i), sex_var, value = T)
  tmp_sex_dat <- codebook %>% filter(toupper(`변수명`) == tmp_sex_var) %>% 
    select(DAT = `테이블명\r\n(영문)`, VAR = `변수명`) 
  tmp_sex_file <- paste0(pheno_path, i, "기/", tmp_sex_dat$DAT, ".csv")
  tmp_sex <- read_csv(tmp_sex_file, show_col_types = FALSE) %>% 
    select(DIST_ID, toupper(tmp_sex_dat$VAR))
  
  tmp_pheno <- tmp_dm %>% left_join(., tmp_sex, by = "DIST_ID")
  
  return(tmp_pheno)
})

pheno_kcdc_2 <- plyr::join_all(extract_dm, by='DIST_ID', type='left') %>%
  select(all_of(colnames(pheno_kcdc_1)))


##### Third Data
file_list <- grep("^\\AS", list.files(pheno_path2), value = TRUE)

pheno_kcdc_3 <- lapply(file_list, function(l) {
  tmp_rdat <- read_delim(paste0(pheno_path2, l), delim= "\t", 
                         show_col_types = FALSE)
  tmp_col <- intersect(colnames(tmp_rdat), c(dm_var, sex_var))
  tmp_dat <- tmp_rdat %>% select(all_of(c("DIST_ID", tmp_col)))
  return(tmp_dat)
})
pheno_kcdc_3 <- pheno_kcdc_3 %>% reduce(left_join, by = "DIST_ID") %>%
  select(all_of(colnames(pheno_kcdc_1)))



##### Age variable (5湲?)
age_dat1 <- read_csv(pheno_kcdc_1_path, show_col_types = FALSE) %>%
  select(DIST_ID = MSKID, AGE = AS5_AGE)

age_dat2 <-  read_csv(paste0(pheno_path, "5기/AS5_01_EXAMINEE.csv"), show_col_types = FALSE) %>% 
  select(DIST_ID, AGE = AS5_AGE)

age_dat3 <- read_delim(paste0(pheno_path2, "AS5_pheno_send_20210825.txt"), delim= "\t", 
                       show_col_types = FALSE) %>%
  select(DIST_ID, AGE = AS5_AGE)


age_dat <- rbind(age_dat1, age_dat2, age_dat3)


### Define DM = TRUE when diagnosed at least once
pheno_kcdc <- rbind(pheno_kcdc_1 %>% mutate(SOURCE = 1), 
                    pheno_kcdc_2 %>% mutate(SOURCE = 2), 
                    pheno_kcdc_3 %>% mutate(SOURCE = 3)) %>% 
  left_join(., age_dat, by = "DIST_ID") %>%
  dplyr::mutate(DM = dplyr::if_any(dm_var, function(x) {x == 2 & !is.na(x)}))

### Check whether SEXs are consitant and return unique value
pheno_kcdc$SEX <- apply(pheno_kcdc %>% select(all_of(sex_var)), 1, 
                        function(x) {
                          re_x <- unlist(ifelse(x %in% c(1, 2, NA), x, NA))
                          unique_sex <- unique(na.omit(re_x))
                          sex <- ifelse(length(unique_sex)!=1, NA, unique_sex)
                          return(sex)
                        })

#pheno_kcdc <- pheno_kcdc %>% mutate(SEX = ifelse(SEX == 1, "M", ifelse(SEX == 2, "F", NA)))

### Age imputation (all the NA values come from source 2)
unique(pheno_kcdc[is.na(pheno_kcdc$AGE), "SOURCE"])

age_na <- pheno_kcdc %>% filter(is.na(AGE))

dat1 <- read_csv(paste0(pheno_path, "1기/AS1_01_EXAMINEE.csv"), show_col_types = FALSE) %>% 
  mutate(ORDER = 1) %>%
  select(DIST_ID, ORDER, AGE = AS1_AGE)

dat2 <- read_csv(paste0(pheno_path, "2기/AS2_01_EXAMINEE.csv"), show_col_types = FALSE) %>% 
  mutate(ORDER = 2) %>%
  select(DIST_ID, ORDER, AGE = AS2_AGE)

dat3 <- read_csv(paste0(pheno_path, "3기/AS3_01_EXAMINEE.csv"), show_col_types = FALSE) %>% 
  mutate(ORDER = 3) %>%
  select(DIST_ID, ORDER, AGE = AS3_AGE)

dat4 <- read_csv(paste0(pheno_path, "4기/AS4_01_EXAMINEE.csv"), show_col_types = FALSE) %>% 
  mutate(ORDER = 4) %>%
  select(DIST_ID, ORDER, AGE = AS4_AGE)

dat <- rbind(dat1, dat2, dat3, dat4)

age_imp <- lapply(age_na$DIST_ID, function(id) {
  dat_age <- dat %>% filter(DIST_ID == id)
  ## Impute by latest age 
  imp_age <- (5-max(dat_age$ORDER))*2 + max(dat_age$AGE)
  res <- data.frame(DIST_ID = id, 
                    AGE_IMP = imp_age)
  
  return(res)
})
age_imp <- do.call("rbind", age_imp)

pheno_kcdc <- pheno_kcdc %>% left_join(., age_imp, by = "DIST_ID") %>%
  mutate(AGE = ifelse(is.na(AGE), AGE_IMP, AGE))


#############################################
##### Bohun
############################################

pheno_bohun <- read_xlsx(paste0(bohun_path, "EPC array total information_rex.xlsx"))
# unique(pheno_bohun$DM) # check whether the data contains all DM subjects


#############################################
##### Combine KCDC and Bohun
############################################
summ_bohun <- pheno_bohun %>% mutate(COHORT = "BOHUN", CASE = 1) %>%
  select(COHORT, CASE, DIST_ID = Sample_ID, SEX = sex, AGE = age, DM)
summ_kcdc <- pheno_kcdc %>% 
  mutate(COHORT = "KCDC", CASE = 0) %>%
  select(COHORT, CASE, DIST_ID, SEX, AGE, DM)

summ_pheno <- rbind(summ_bohun, summ_kcdc)

detach(package:plyr) 
summ_pheno %>% 
  dplyr::group_by(COHORT, CASE, SEX) %>% 
  summarize(N = length(unique(DIST_ID)))

write.table(summ_pheno, paste0(out_path, "01_01_DataPrep_Pheno_age.csv"), 
            sep = ",",
            row.names = F, col.names = T)




