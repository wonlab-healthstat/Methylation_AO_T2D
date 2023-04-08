####################################
### Study : Bohun_Methylation
### Purpose of program : Convert vcf files to bed/bim/fam files
### Pre-required program : NA
### Programmed by : Sujin Seo
### Output : PLINK format file
### Draft Date : 2022-04-12
### Revision : NA
### Software (version) : R (v4.1.1)
#####################################

#####################################
########## Load Required Libraries
#####################################
library(dplyr)
library(parallel)
library(stringr)

#####################################
########## Set Parameters
#####################################
bh_dm_dir <- "/data/bohun_methylation/4171/질병관리본부데이터/KCHIP/"
dr_dir <- "/home2/sjseo/Bohun/SNParray/data/"

ncores <- 20

#####################################
########## Data Preparation
#####################################

mclapply(1:22, function(i) {
  
  print(paste0("Chr", i, " start!!!\n"))
  
  ### Filter (Imputation quality)
  filter_cmd1 <- paste("vcftools --gzvcf", paste0(bh_dm_dir, "BH_chr", i, "_20200813.vcf.gz"),
                       "--remove-filtered-all --recode-INFO-all --recode",
                       "--out", paste0(dr_dir, "BH_DM/InfoFilter/BH_chr", i))
  system(filter_cmd1)
 
  ### Convert to bed/bim/fam
  cvt_plk_cmd <- paste("plink --vcf", paste0(dr_dir, "BH_DM/InfoFilter/BH_chr", i, ".recode.vcf"), 
                       "--make-bed --out", paste0(dr_dir, "BH_DM/Chr", i))
  system(cvt_plk_cmd)
}, mc.cores = ncores)
### InfoFilter folder was deleted after the running of 01_01_DataPrep_DM.R code.

### Merge to one file
file_list <- c()
for(i in 1:22) {
  file_list <- c(file_list, paste0(paste0(dr_dir, "BH_DM/Chr", i), c(".bed", ".bim", ".fam"), collapse = " "))
}
write.table(file_list, paste0(dr_dir, "BH_DM/file_list.txt"), sep = " ", col.names = F, row.names = F, quote = F)
merge_cmd <- paste("plink --merge-list", paste0(dr_dir, "BH_DM/file_list.txt"), "--out", paste0(dr_dir, "BH_DM/BH_DM_merged"))
system(merge_cmd)

rm_cmd1 <- paste("rm", paste0(dr_dir, "BH_DM/Chr*"))
rm_cmd2 <- paste("rm", paste0(dr_dir, "BH_DM/file_list.txt"))
system(rm_cmd1); systerm(rm_cmd2)

### ID correction
fam <- read.table(paste0(dr_dir, "BH_DM/BH_DM_merged.fam"), header = F, stringsAsFactors = F)
new_id <- paste0(fam$V1, "_", str_pad(fam$V2, 3, pad = "0"))
fam$V1 <- fam$V2 <- new_id
write.table(fam, paste0(dr_dir, "BH_DM/BH_DM_merged.fam"), sep = "\t", col.names = F, row.names = F, quote = F)


