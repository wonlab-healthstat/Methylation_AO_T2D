#####################################
### Study : Bohun-Methylation
### Purpose of program : Create pheno data(DM, SEX)
### Pre-required program : NA
### Output files : KoGES_DM_SEX.csv
### Programmed by : Sujin Seo
### Draft Date : 2022-04-25
### Revision : NA
### Software (version) : R (v3.6.3)
#####################################

#####################################
##### Load Required Libraries
#####################################
library(readxl)
library(dplyr)
library(tidyverse)

func_comb_dm <- function(dat, idx) {
  comb_dat <- apply(dat, 1, function(x) {
    na_idx <- which(is.na(x[idx]) | x[idx] %in% c(66666, 77777, 99999))
    if(length(na_idx) > 0) {
      inc_idx <- idx[-na_idx]
    } else {
      inc_idx <- idx
    }
    val <- max(as.numeric(x[inc_idx]), na.rm = T)
    ifelse(val %in% c(1, 2), val, NA)
  })
  return(comb_dat)
}

func_comb_sex <- function(dat, idx) {
  comb_dat <- apply(dat, 1, function(x) {
    na_idx <- which(is.na(x[idx]) | x[idx] %in% c(66666, 77777, 99999))
    if(length(na_idx) > 0) {
      inc_idx <- idx[-na_idx]
    } else {
      inc_idx <- idx
    }
    uniq_sex <- unique(as.numeric(x[inc_idx]))
    val <- ifelse(length(uniq_sex) == 1, uniq_sex, NA)
    ifelse(val %in% c(1, 2), val, NA)
  })
  return(comb_dat)
}


#####################################
##### Set parameters
#####################################

##### Directories of KARE codingbook
kare_dir <- "/data2/KoGES_2021/KoGES_AS_KARE/Phenotype/CodingBook/"

##### Directories of NC codingbook
nc_dir <- "/data2/KoGES_2021/KoGES_NC_CAVAS/Phenotype/CodingBook/"

##### Directories of CT codingbook
ct_dir <- "/data2/KoGES_2021/KoGES_CT_HEXA/Phenotype/CodingBook/"

##### Directories of bohun
bh_dir <- "/data/bohun_methylation/"

##### Derived directory
dr_dir <- "/home2/sjseo/Bohun/SNParray/derived_data/"

###########################################
##### Kare Data
###########################################

### Diabetes variable name for each cohort
dm_var <- c("AS1_PdDm", "AS2_PdDm", "AS3_PdFDm", "AS4_DMDIAG", 
            "AS5_DMDIAG", "AS6_DMDIAG", "AS7_DMDIAG", "AS8_DMDIAG")

### Sex variable name for each cohort
sex_var <- paste0("AS", 1:8, "_SEX")

### List of codebooks
file_lst <- list.files(kare_dir)

### Extract DM/SEX
kare_pheno <- lapply(1:8, function(i) {
  file_nm <- grep(paste0("^(AS", i, "_)"), file_lst, value = T)
  codebook <- read_excel(paste0(kare_dir, file_nm),
                         sheet = "공개코드북", skip = 2)
  
  # DM Information
  tmp_dm_var <- grep(paste0("^\\AS", i), dm_var, value = T)
  tmp_dm_dat <- codebook %>% filter(`변수명` == tmp_dm_var) %>% 
    dplyr::select(DAT = `테이블명\r\n(영문)`, VAR = `변수명`) 
  
   
  # Sex Information
  tmp_sex_var <- grep(paste0("^\\AS", i), sex_var, value = T)
  tmp_sex_dat <- codebook %>% filter(toupper(`변수명`) == tmp_sex_var) %>% 
    dplyr::select(DAT = `테이블명\r\n(영문)`, VAR = `변수명`) 
  
  
  # File list of phenotype data
  tmp_dir <- grep(paste0("^(AS", i, "_)"), list.files(paste0(kare_dir, "../")), value = T)
  
  tmp_dm_file <- paste0(kare_dir, "../", tmp_dir, "/", tmp_dm_dat$DAT, ".csv")
  tmp_dm <- read_csv(tmp_dm_file, show_col_types = FALSE) %>% 
    dplyr::select(DIST_ID, DM = toupper(tmp_dm_dat$VAR))
  
  tmp_sex_file <- paste0(kare_dir, "../", tmp_dir, "/", tmp_sex_dat$DAT, ".csv")
  tmp_sex <- read_csv(tmp_sex_file, show_col_types = FALSE) %>% 
    dplyr::select(DIST_ID, SEX = toupper(tmp_sex_dat$VAR))
  
  tmp_pheno <- tmp_dm %>% left_join(., tmp_sex, by = "DIST_ID")
  
  return(tmp_pheno)
})


### Combine results
kare_pheno_comb <- Reduce(function(x, y) full_join(x, y, by = "DIST_ID"), kare_pheno) 

dm_idx <- grep("^(DM)", colnames(kare_pheno_comb))
sex_idx <- grep("^(SEX)", colnames(kare_pheno_comb))

kare_pheno_comb$DM <- func_comb_dm(dat = kare_pheno_comb, idx = dm_idx)
kare_pheno_comb$SEX <- func_comb_sex(dat = kare_pheno_comb, idx = sex_idx)
kare_pheno_comb <- kare_pheno_comb[, c("DIST_ID", "DM", "SEX")]


###########################################
##### NC Data
###########################################

### Diabetes variable name for each cohort
dm_var <- paste0("NC", 1:5, "_DM")

### Sex variable name for each cohort
sex_var <- paste0("NC", 1:5, "_SEX")

### List of codebooks
file_lst <- list.files(nc_dir)

### Extract DM/SEX
nc_pheno <- lapply(1:5, function(i) {
  file_nm <- grep(paste0("^(NC", i, "_)"), file_lst, value = T)
  codebook <- read_excel(paste0(nc_dir, file_nm),
                         sheet = 2, skip = 2)
  
  # DM Information
  tmp_dm_var <- grep(paste0("^\\NC", i), dm_var, value = T)
  tmp_dm_dat <- codebook %>% filter(`변수명` == tmp_dm_var) %>% 
    dplyr::select(DAT = `테이블명\r\n(영문)`, VAR = `변수명`) 
 
  # Sex Information
  tmp_sex_var <- grep(paste0("^\\NC", i), sex_var, value = T)
  tmp_sex_dat <- codebook %>% filter(toupper(`변수명`) == tmp_sex_var) %>% 
    dplyr::select(DAT = `테이블명\r\n(영문)`, VAR = `변수명`) 
  
 
  # File list of phenotype data
  tmp_dir <- grep(paste0("^(NC", i, "_)"), list.files(paste0(nc_dir, "../")), value = T)
  
  tmp_dm_file <- paste0(nc_dir, "../", tmp_dir, "/", tmp_dm_dat$DAT, ".csv")
  tmp_dm <- read_csv(tmp_dm_file, show_col_types = FALSE) %>% 
    dplyr::select(DIST_ID, DM = toupper(tmp_dm_dat$VAR))
  
  tmp_sex_file <- paste0(nc_dir, "../", tmp_dir, "/", tmp_sex_dat$DAT, ".csv")
  tmp_sex <- read_csv(tmp_sex_file, show_col_types = FALSE) %>% 
    dplyr::select(DIST_ID, SEX = toupper(tmp_sex_dat$VAR))
  
  tmp_pheno <- tmp_dm %>% left_join(., tmp_sex, by = "DIST_ID")
  
  return(tmp_pheno)
})

### Combine results
nc_pheno_comb <- Reduce(function(x, y) full_join(x, y, by = "DIST_ID"), nc_pheno) 

dm_idx <- grep("^(DM)", colnames(nc_pheno_comb))
sex_idx <- grep("^(SEX)", colnames(nc_pheno_comb))

nc_pheno_comb$DM <- func_comb_dm(dat = nc_pheno_comb, idx = dm_idx)
nc_pheno_comb$SEX <- func_comb_sex(dat = nc_pheno_comb, idx = sex_idx)
nc_pheno_comb <- nc_pheno_comb[, c("DIST_ID", "DM", "SEX")]


###########################################
##### CT Data
###########################################

### Diabetes variable name for each cohort
dm_var <- paste0("CT", 1:2, "_DM")

### Sex variable name for each cohort
sex_var <- paste0("CT", 1:2, "_SEX")

### List of codebooks
file_lst <- list.files(ct_dir)

### Extract DM/SEX
ct_pheno <- lapply(1:2, function(i) {
  file_nm <- grep(paste0("^(CT", i, "_)"), file_lst, value = T)
  codebook <- read_excel(paste0(ct_dir, file_nm),
                         sheet = 2, skip = 2)
  
  # DM Information
  tmp_dm_var <- grep(paste0("^\\CT", i), dm_var, value = T)
  tmp_dm_dat <- codebook %>% filter(`변수명` == tmp_dm_var) %>% 
    dplyr::select(DAT = `테이블명\r\n(영문)`, VAR = `변수명`) 
  
  # Sex Information
  tmp_sex_var <- grep(paste0("^\\CT", i), sex_var, value = T)
  tmp_sex_dat <- codebook %>% filter(toupper(`변수명`) == tmp_sex_var) %>% 
    dplyr::select(DAT = `테이블명\r\n(영문)`, VAR = `변수명`) 
  
  # File list of phenotype data
  tmp_dir <- grep(paste0("^(CT", i, "_)"), list.files(paste0(ct_dir, "../")), value = T)
  
  tmp_dm_file <- paste0(ct_dir, "../", tmp_dir, "/", tmp_dm_dat$DAT, ".csv")
  tmp_dm <- read_csv(tmp_dm_file, show_col_types = FALSE) %>% 
    dplyr::select(DIST_ID, DM = toupper(tmp_dm_dat$VAR))
  
  tmp_sex_file <- paste0(ct_dir, "../", tmp_dir, "/", tmp_sex_dat$DAT, ".csv")
  tmp_sex <- read_csv(tmp_sex_file, show_col_types = FALSE) %>% 
    dplyr::select(DIST_ID, SEX = toupper(tmp_sex_dat$VAR))
  
  tmp_pheno <- tmp_dm %>% left_join(., tmp_sex, by = "DIST_ID")
  
  return(tmp_pheno)
})


### Combine results
ct_pheno_comb <- Reduce(function(x, y) full_join(x, y, by = "DIST_ID"), ct_pheno) 

dm_idx <- grep("^(DM)", colnames(ct_pheno_comb))
sex_idx <- grep("^(SEX)", colnames(ct_pheno_comb))

ct_pheno_comb$DM <- func_comb_dm(dat = ct_pheno_comb, idx = dm_idx)
ct_pheno_comb$SEX <- func_comb_sex(dat = ct_pheno_comb, idx = sex_idx)
ct_pheno_comb <- ct_pheno_comb[, c("DIST_ID", "DM", "SEX")]



###########################################
##### Bohun Data
###########################################
##### Vitamin D
bh_vitamin <- read_xlsx(paste0(bh_dir, "bohun_vitd_final.xlsx")) %>% #503
  mutate(DM = sapply(MEDICAL_HISTORY, function(x) ifelse(grepl("(당뇨)", x), 2, 1))) %>%
  select(DIST_ID = DONOR_ID, DM, SEX = "성별")

##### Muscle 
bh_muscle <- read_xlsx(paste0(bh_dir, "bohun_muscle_study_clinic_data.xlsx")) %>% #1,695
  mutate(DM = sapply(MEDICAL_HISTORY, function(x) ifelse(grepl("(당뇨)", x), 2, 1))) %>%
  select(DIST_ID = DONOR_ID, DM, SEX = "성별")

###### Merge above to datasets
bh_pheno <- unique(rbind(bh_vitamin, bh_muscle)) %>%
  mutate(SEX =as.numeric(ifelse(SEX == "남자", 1, ifelse(SEX == "여자", 2, SEX))))



###########################################
##### Combine KARE, NC, CT, Bohun
###########################################
nrow(kare_pheno_comb) + nrow(ct_pheno_comb) + nrow(nc_pheno_comb)
length(unique(c(kare_pheno_comb$DIST_ID, nc_pheno_comb$DIST_ID, ct_pheno_comb$DIST_ID)))

comb_dat <- rbind(kare_pheno_comb, nc_pheno_comb, ct_pheno_comb, bh_pheno) #81,838
write.table(comb_dat, 
            file.path(dr_dir, "02_02_Pheno_DM_SEX.csv"),
            sep = ",", col.names = T, row.names = F, quote = F)
# comb_dat <- read.table(file.path(dr_dir, "02_02_Pheno_DM_SEX.csv"), sep = ",", header = T)


###########################################
##### Test data should be excluded
##### Test data is subset of whole data
###########################################

test_pheno <- read.table(file.path(dr_dir, "02_01_IDmatched.fam"), sep = ' ', header = F, stringsAsFactors = F)
match_id <- read.table(file.path(dr_dir, "02_01_IDmatch/Matched_ID.txt"), sep = "\t", header = T, stringsAsFactors = F)
test_id <- match_id[which(match_id$Geno_Bohun %in% test_pheno$V2), "Geno_Wonlab"]
ref_id <- comb_dat %>% filter(!(DIST_ID %in% test_id)) %>%
  select(FID = DIST_ID, IID = DIST_ID, DM, SEX)

write.table(ref_id, 
            file.path(dr_dir, "02_02_Pheno_DM_SEX_ref.txt"),
            sep = " ", col.names = T, row.names = F, quote = F, na = "-9")

