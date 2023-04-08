#####################################
### Study : Bohun-Methylation
### Purpose of program : Baseline demographics
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
library(readxl)
library(dplyr)
library(readr)
library(plyr)
library(purrr)


#####################################
##### Set Working Directories
#####################################

dr_dir <- "/home2/sjseo/Bohun/Methylation/derived_data/"


#############################################
##### Target variables
############################################

##### Only refer to 5th cohort
var1 <- c("AGE", "SEX", "BMI", "EDATE5", "HBA1C", "CREATININE")

###### Refer 5th and previous cohort
#como <- c("MI", "HTN", "LCA", "GCA", "HCC", "COLCA", "PACA", "UTCA", "BRCA", "CEVA")
como <- c("MI", "HTN", "LCA", "GCA", "HCC", "COLCA", "PACA", "CEVA")
como_crnt <- paste0(como, "CU")
var2 <- c("DMDIAG", "DMDIAGYR", "DMDIAGMO", como, como_crnt)

###### All the targeted variables
var <- c(var1, var2)
missing <- c(66666, 77777, 99999)


#############################################
##### KCDC 4171
############################################

##### Cohort 1 (4171)
coh1 <- read_csv("/data/bohun_methylation/4171/질병관리본부데이터/methlaytion_dat.csv", show_col_types = FALSE)
colnames(coh1)[which(colnames(coh1) == "AS5_CREATININ")] <- "AS5_CREATININE"

##### Check the existence of target variables
sapply(var, function(v) grep(paste0("_(", v, ")$"), colnames(coh1), value = T))
tmp_var2 <- as.vector(unlist(sapply(var2, function(v) grep(paste0("_(", v, ")$"), colnames(coh2), value = T))))


##### Common variables
coh1_tmp1 <- coh1 %>% select(all_of(c("MSKID", paste0("AS5_", var1)))) %>%
  mutate(tmp_eGRF = 175 * (AS5_CREATININE^-1.154) * (AS5_AGE^-0.203),
         AS5_eGRF = ifelse(AS5_SEX == 2, tmp_eGRF*0.742, tmp_eGRF),
         BIRTH = as.Date(paste0(as.numeric(substr(AS5_EDATE5, 1, 4)) - AS5_AGE, 
                                "-", substr(AS5_EDATE5, 5, 6),
                                "-01"), format = "%Y-%m-%d")) %>%
  select(-tmp_eGRF)

#### Diabetes
coh1_tmp2 <- coh1 %>% select(all_of(c("MSKID", tmp_var2))) %>%
  left_join(., coh1_tmp1 %>% select(MSKID, AS5_EDATE5, BIRTH), by = "MSKID") %>% 
  mutate(DM = ifelse(AS4_DMDIAG == 2 | AS5_DMDIAG == 2, 2, 1),
         DMDIAGDATE4 = as.Date(ifelse(AS4_DMDIAGYR %in% missing | AS4_DMDIAGMO %in% missing, NA, 
                                      paste0(AS4_DMDIAGYR,"-", str_pad(AS4_DMDIAGMO, 2, pad = "0"), "-01")),
                               format = "%Y-%m-%d"),
         DMDIAGDATE5 = as.Date(ifelse(AS5_DMDIAGYR %in% missing | AS5_DMDIAGMO %in% missing, NA, 
                                      paste0(AS5_DMDIAGYR,"-", str_pad(AS5_DMDIAGMO, 2, pad = "0"), "-01")),
                               format = "%Y-%m-%d"),
         DMDIAGDATE = pmin(DMDIAGDATE4, DMDIAGDATE5, na.rm = T),
         DMONSET = as.numeric(difftime(DMDIAGDATE, BIRTH, unit = "weeks")/52.25),
         AS5_EDATE5 = as.Date(paste0(as.numeric(substr(AS5_EDATE5, 1, 4)), 
                                     "-", substr(AS5_EDATE5, 5, 6),
                                     "-01"), format = "%Y-%m-%d"),
         DMDUR = as.numeric(difftime(AS5_EDATE5, DMDIAGDATE, unit = "weeks")/52.25)) %>%
  select(MSKID, DM, DMDIAGDATE, DMONSET, DMDUR)


##### Comobidity
coh1_tmp3 <- lapply(como, function(v) {
  
  print(v)
  
  vcu <- paste0(v, "CU")
  var2 <- rep(c(v, vcu), 2)
  order <- rep(c("AS4", "AS5"), each = 2)
  target_var <- paste(order , var2, sep = "_")
  target_var <- c("MSKID", intersect(target_var, tmp_var2))

  if(paste0("AS5_", vcu) %in% target_var) {
    
    tmp_como <- data.frame(coh1) %>% select(all_of(target_var)) %>%
      mutate(tmp_as4_v = ifelse(get(paste0("AS4_", v)) %in% missing , NA,  get(paste0("AS4_", v))),
             tmp_as5_v = ifelse(get(paste0("AS5_", v)) %in% missing , NA,  get(paste0("AS5_", v))),
             tmp_as4_vcu = ifelse(get(paste0("AS4_", vcu)) %in% missing , NA,  get(paste0("AS4_", vcu))),
             tmp_as5_vcu = ifelse(get(paste0("AS5_", vcu)) %in% missing , NA,  get(paste0("AS5_", vcu))))
    
    tmp_como[,v] <- with(tmp_como, pmax(tmp_as4_v, tmp_as5_v, na.rm = T))
    tmp_como[,vcu] <- with(tmp_como, pmin(tmp_as4_vcu, tmp_as5_vcu, na.rm = T))         
    tmp_como[,vcu] <- with(tmp_como, ifelse(get(vcu) == 1, 1, 2))         
    tmp_como <- tmp_como[, c("MSKID", v, vcu)] 

  } else {
    
    tmp_como <- data.frame(coh1) %>% select(all_of(target_var)) %>%
      mutate(tmp_as4_v = ifelse(get(paste0("AS4_", v)) %in% missing , NA,  get(paste0("AS4_", v))),
             tmp_as5_v = ifelse(get(paste0("AS5_", v)) %in% missing , NA,  get(paste0("AS5_", v))),
             tmp_as4_vcu = ifelse(get(paste0("AS4_", vcu)) %in% missing , NA,  get(paste0("AS4_", vcu))))
    
    tmp_como[,v] <- with(tmp_como, pmax(tmp_as4_v, tmp_as5_v, na.rm = T))
    tmp_como[,vcu] <- tmp_como$tmp_as4_vcu         
    tmp_como[,vcu] <- with(tmp_como, ifelse(get(vcu) == 1, 1, 2))         
    tmp_como <- tmp_como[, c("MSKID", v, vcu)] 
    
  }
  
  return(tmp_como)
})
coh1_tmp3 <- plyr::join_all(coh1_tmp3, by='MSKID') 

coh1_res <- plyr::join_all(list(coh1_tmp1, coh1_tmp2, coh1_tmp3), by='MSKID') 
colnames(coh1_res)[1] <- "DIST_ID"


#############################################
##### KCDC 5566
############################################

##### Cohort 2 (5566)
coh2 <- lapply(1:5, function(i) {
  tmp <- read.table(paste0("/data/bohun_methylation/5566/메틸레이션_질본_정상_추가/AS", i, "_pheno_send_20210825.txt"), 
                    sep = "\t", header = T, stringsAsFactors = F)
})
coh2 <- plyr::join_all(coh2, by='DIST_ID')
colnames(coh2)[which(colnames(coh2) == "AS5_CREATININ")] <- "AS5_CREATININE"

##### Check the existence of target variables
sapply(var, function(v) grep(paste0("_(", v, ")$"), colnames(coh2), value = T))
tmp_var2 <- unlist(sapply(var2, function(v) grep(paste0("_(", v, ")$"), colnames(coh2), value = T)))

##### Common variables
coh2_tmp1 <- coh2 %>% select(all_of(c("DIST_ID", paste0("AS5_", var1)))) %>%
  mutate(tmp_eGRF = 175 * (AS5_CREATININE^-1.154) * (AS5_AGE^-0.203),
         AS5_eGRF = ifelse(AS5_SEX == 2, tmp_eGRF*0.742, tmp_eGRF),
         BIRTH = as.Date(paste0(as.numeric(substr(AS5_EDATE5, 1, 4)) - AS5_AGE, 
                                "-", substr(AS5_EDATE5, 5, 6),
                                "-01"), format = "%Y-%m-%d")) %>%
  select(-tmp_eGRF)


#### Diabetes
coh2_tmp2 <- coh2 %>% select(all_of(c("DIST_ID", tmp_var2))) %>%
  left_join(., coh2_tmp1 %>% select(DIST_ID, AS5_EDATE5, BIRTH), by = "DIST_ID") %>% 
  mutate(DM = ifelse(AS4_DMDIAG == 2 | AS5_DMDIAG == 2, 2, 1),
         DMDIAGDATE4 = as.Date(ifelse(AS4_DMDIAGYR %in% missing | AS4_DMDIAGMO %in% missing, NA, 
                                      paste0(AS4_DMDIAGYR,"-", str_pad(AS4_DMDIAGMO, 2, pad = "0"), "-01")),
                               format = "%Y-%m-%d"),
         DMDIAGDATE5 = as.Date(ifelse(AS5_DMDIAGYR %in% missing | AS5_DMDIAGMO %in% missing, NA, 
                                      paste0(AS5_DMDIAGYR,"-", str_pad(AS5_DMDIAGMO, 2, pad = "0"), "-01")),
                               format = "%Y-%m-%d"),
         DMDIAGDATE = pmin(DMDIAGDATE4, DMDIAGDATE5, na.rm = T),
         DMONSET = as.numeric(difftime(DMDIAGDATE, BIRTH, unit = "weeks")/52.25),
         AS5_EDATE5 = as.Date(paste0(as.numeric(substr(AS5_EDATE5, 1, 4)), 
                                     "-", substr(AS5_EDATE5, 5, 6),
                                     "-01"), format = "%Y-%m-%d"),
         DMDUR = as.numeric(difftime(AS5_EDATE5, DMDIAGDATE, unit = "weeks")/52.25)) %>%
  select(DIST_ID, DM, DMDIAGDATE, DMONSET, DMDUR)


##### Comobidity
coh2_tmp3 <- lapply(como, function(v) {
  
  print(v)
  
  vcu <- paste0(v, "CU")
  var2 <- rep(c(v, vcu), 2)
  order <- rep(c("AS4", "AS5"), each = 2)
  target_var <- paste(order , var2, sep = "_")
  target_var <- c("DIST_ID", intersect(target_var, tmp_var2))
  
  if(paste0("AS5_", vcu) %in% target_var) {
    
    tmp_como <- data.frame(coh2) %>% select(all_of(target_var)) %>%
      mutate(tmp_as4_v = ifelse(get(paste0("AS4_", v)) %in% missing , NA,  get(paste0("AS4_", v))),
             tmp_as5_v = ifelse(get(paste0("AS5_", v)) %in% missing , NA,  get(paste0("AS5_", v))),
             tmp_as4_vcu = ifelse(get(paste0("AS4_", vcu)) %in% missing , NA,  get(paste0("AS4_", vcu))),
             tmp_as5_vcu = ifelse(get(paste0("AS5_", vcu)) %in% missing , NA,  get(paste0("AS5_", vcu))))
    
    tmp_como[,v] <- with(tmp_como, pmax(tmp_as4_v, tmp_as5_v, na.rm = T))
    tmp_como[,vcu] <- with(tmp_como, pmin(tmp_as4_vcu, tmp_as5_vcu, na.rm = T))         
    tmp_como[,vcu] <- with(tmp_como, ifelse(get(vcu) == 1, 1, 2))         
    tmp_como <- tmp_como[, c("DIST_ID", v, vcu)] 
    
  } else {
    
    tmp_como <- data.frame(coh2) %>% select(all_of(target_var)) %>%
      mutate(tmp_as4_v = ifelse(get(paste0("AS4_", v)) %in% missing , NA,  get(paste0("AS4_", v))),
             tmp_as5_v = ifelse(get(paste0("AS5_", v)) %in% missing , NA,  get(paste0("AS5_", v))),
             tmp_as4_vcu = ifelse(get(paste0("AS4_", vcu)) %in% missing , NA,  get(paste0("AS4_", vcu))))
    
    tmp_como[,v] <- with(tmp_como, pmax(tmp_as4_v, tmp_as5_v, na.rm = T))
    tmp_como[,vcu] <- tmp_como$tmp_as4_vcu         
    tmp_como[,vcu] <- with(tmp_como, ifelse(get(vcu) == 1, 1, 2))         
    tmp_como <- tmp_como[, c("DIST_ID", v, vcu)] 
    
  }
  
  return(tmp_como)
})
coh2_tmp3 <- plyr::join_all(coh2_tmp3, by='DIST_ID') 

coh2_res <- plyr::join_all(list(coh2_tmp1, coh2_tmp2, coh2_tmp3), by='DIST_ID') 


#############################################
##### Combine 4171/5566 cohorts and export 
############################################

res_dat <- rbind(coh1_res, coh2_res)
write.table(res_dat, paste0(dr_dir, "01_02_Baseline_KARE.csv"), sep = ",", col.names = T, row.names = F, quote = F)
