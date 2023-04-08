#####################################
### Study : Bohun-Methylation
### Purpose of prgram : Read *idat and convert to R object
### Pre-required program : NA
### Output files : Files in "meth_rdata" folder
### Programmed by : Sujin Seo
### Draft Date : 2021-09-30
### Revision : NA
### Software (version) : R (v3.6.3)
#####################################


#####################################
##### Load Required Libraries
#####################################
library(ewastools)
library(dplyr)
library(readr)
library(parallel)


#####################################
##### Set Parameters
#####################################
##### *.idat folder path
path_bh1_1 <- "/data/bohun_methylation/4171/multiomics/멀티오믹스를 이용한 근감소증 바이오마커발굴/5. 고엽제 당뇨 EWAS 비고엽제당뇨 EWAS Epic 850K/1차 당뇨 메틸레이션 53명/EPIC array 1차 53명/191004-190924AM-01_Bohun_SJH_EPIC_48_raw/(1)Raw.data/idat/"
path_bh1_2 <- "/data/bohun_methylation/4171/multiomics/멀티오믹스를 이용한 근감소증 바이오마커발굴/5. 고엽제 당뇨 EWAS 비고엽제당뇨 EWAS Epic 850K/1차 당뇨 메틸레이션 53명/EPIC array 1차 53명/191030-191023AM-01_Bohun_SJH_EPIC_5_raw/(1)Raw.data/idat/"
path_bh2 <- "/data/bohun_methylation/4171/multiomics/멀티오믹스를 이용한 근감소증 바이오마커발굴/5. 고엽제 당뇨 EWAS 비고엽제당뇨 EWAS Epic 850K/2차 당뇨 메틸레이션 72명/idat/"
path_kcdc1 <- "/data/bohun_methylation/4171/질병관리본부데이터/Methylation"
path_kcdc2 <- "/data/bohun_methylation/5566/메틸레이션_질본_정상_추가/Methylation_850K_863"
path <- c(path_bh1_1, path_bh1_2, path_bh2, path_kcdc1, path_kcdc2)


##### Derived Data Folder
dr_dir <- "/home2/sjseo/Bohun/Methylation/derived_data/"
dir.create(file.path(dr_dir,"meth_rdata"))

##### Number of Cores
ncore <- 4


#####################################
##### Read Data & Data Preparation
#####################################

##### Pheno Data
exp_var <- read.table(paste0(dr_dir, "DataPrep_Pheno.csv"), sep = ",",
                      header = T, stringsAsFactors=FALSE)


##### Sample 
##### IDAT of bohun has sentrix id as a sample_id, so combine two information
bh_smp_file <- c("/data/bohun_methylation/4171/multiomics/멀티오믹스를 이용한 근감소증 바이오마커발굴/5. 고엽제 당뇨 EWAS 비고엽제당뇨 EWAS Epic 850K/1차 당뇨 메틸레이션 53명/EPIC array 1차 53명/Sample.Table 1.txt",
                 "/data/bohun_methylation/4171/multiomics/멀티오믹스를 이용한 근감소증 바이오마커발굴/5. 고엽제 당뇨 EWAS 비고엽제당뇨 EWAS Epic 850K/1차 당뇨 메틸레이션 53명/EPIC array 1차 53명/Sample.Table 2.txt",
                 "/data/bohun_methylation/4171/multiomics/멀티오믹스를 이용한 근감소증 바이오마커발굴/5. 고엽제 당뇨 EWAS 비고엽제당뇨 EWAS Epic 850K/2차 당뇨 메틸레이션 72명/Sample.Table.txt")

bh_smp_file <- do.call("rbind",
                       lapply(bh_smp_file, function(x) dat <- read_delim(x, delim = "\t", show_col_types = FALSE)))
bh_smp_file <- bh_smp_file %>%
  mutate(smp = paste0(`Sentrix Barcode`, "_", `Sample Section`)) %>%
  select(`Sample ID`, smp)

exp_var <- exp_var %>% 
  left_join(., bh_smp_file, by = c("DIST_ID" = "Sample ID")) %>%
  mutate(DIST_ID, Sample_ID = ifelse(is.na(smp), DIST_ID, smp)) %>%
  select(-c(DIST_ID, smp))

exp_var %>% 
  dplyr::group_by(COHORT, CASE, DM, SEX) %>% 
  summarize(N = length(unique(Sample_ID)))


write.table(exp_var, paste0(dr_dir, "02_01_Read_idat_Pheno_Anal.csv"), sep = ",", col.names = T, row.names = F)

##### IDAT Data
slides <- lapply(path, function(p) {
  file_list <- list.files(p, recursive=TRUE)
  slide <- unique(sapply(file_list, 
                         function(f) strsplit(f, "(_Red*)|(_Grn*)")[[1]][1]))
  
  return(data.frame(smp = slide, 
                    path = paste0(p, slide),
                    stringsAsFactors = F))
  
})
slides <- do.call("rbind", slides) %>% filter(smp %in% exp_var$Sample_ID)

read_export_idat <- function(i) {
  s <- slides[i, "smp"]
  p <- slides[i, "path"]
  tmp_dat <- read_idats(p, quiet = TRUE)
  save(tmp_dat, file = paste0(dr_dir, "/meth_rdata/", s, ".RData"))
  rm(list = c("s", "p", "tmp_dat"))
}
mclapply(1:nrow(slides), read_export_idat, mc.cores = ncore)
# Delete meth_rdata folder after run 03_01_QC.R for memory saving

