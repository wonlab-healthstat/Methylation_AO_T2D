#####################################
### Study : Bohun-Methylation
### Purpose of program : GWAS with discovery data
### Pre-required program : NA
### Output files : 
### Programmed by : Sujin Seo
### Draft Date : 2022-05-02
### Revision : NA
### Software (version) : R (v3.6.3)
#####################################

#####################################
##### Load Required Libraries
#####################################
library(parallel)
library(qqman)
library(dplyr)

#####################################
##### Set parameters
#####################################

##### Derived directory
dr_dir <- "/home2/sjseo/Bohun/SNParray/derived_data/03_02_Logit_DM/"

##### Number of cores
ncores <- detectCores()


#####################################
##### Logistic regression (PLINK)
#####################################

mc_cmd <- mclapply(1:22, function(chr) {
  gwas_cmd <- paste("plink --bfile", paste0("/data4/ardo/Data2/20201102_ImpuatationTest/NC_CT_KARE_Bohun_Merge/Merge_kareAffy_ncAffy_ctAffy_bhKchip_kncKchip/4.snpQC_genomafhwe/",
                                            "cleaned_chr", chr, "_Merge_kareAffy_ncAffy_ctAffy_bhKchip_kncKchip"),
                    "--keep", paste0(dr_dir, "../03_01_Pheno_DM_SEX_ref.txt"),
                    "--covar", paste0(dr_dir, "../03_01_Pheno_DM_SEX_ref.txt"), 
                    "--covar-name SEX", 
                    "--pheno", paste0(dr_dir, "../03_01_Pheno_DM_SEX_ref.txt"),
                    "--pheno-name DM",
                    "--logistic",
                    "--out", paste0(dr_dir, "Logit_DM_chr", chr))
  system(gwas_cmd)
}, mc.cores = ncores)


#####################################
##### Combine the above result
#####################################

logit_res <- do.call(rbind, 
                     mclapply(1:22, function(chr) 
                       read.table(paste0(dr_dir, "Logit_DM_chr", chr, ".assoc.logistic"), header = T, stringsAsFactors = F),
                     mc.cores = ncores)) %>%
  filter(TEST == "ADD") %>%
  mutate(BETA = log(OR), SE = BETA/STAT)
  
#logit_res <- read.table(paste0(dr_dir, "Logit_DM.assoc.logistic"), header = T, stringsAsFactors = F) %>%
#  filter(TEST == "ADD") %>%
#  mutate(BETA = log(OR)) 

#### Reference allele
logit_res$A2 <- unlist(mclapply(1:nrow(logit_res), function(i) {
  x <- logit_res[i, ]
  a2 <- setdiff(strsplit(x$SNP, ":")[[1]][3:4], x$A1)
}, mc.cores = ncores))

logit_res <- rename(logit_res, REF = A2, ALT = A1, N = NMISS)
write.table(logit_res, paste0(dr_dir, "../03_03_PRS_Wonlab/DM_discovery.sumstats"), sep = "\t", col.names = T, row.names = F, quote = F)
#logit_res <- read.table(paste0(dr_dir, "../03_03_PRS_Wonlab/DM_discovery.sumstats"), sep = "\t", header = T, stringsAsFactors = F)

##### Manhattan Plot
manhattan(logit_res, chr = "CHR", bp = "BP", p = "P", snp = "SNP",
  col = c("gray10", "gray60"),
  suggestiveline = -log10(1e-05),
  genomewideline = -log10(5e-08))

##### QQ Plot
qq(logit_res$P)


