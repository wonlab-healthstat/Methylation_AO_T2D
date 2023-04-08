####################################
### Study : Bohun_Methylation
### Purpose of program : Combine DM and Normal data
### Pre-required program : NA
### Programmed by : Sujin Seo
### Output : PLINK format file
### Draft Date : 2022-02-25
### Revision : NA
### Software (version) : R (v4.1.1)
#####################################

#####################################
########## Load Required Libraries
#####################################
library(dplyr)
library(parallel)


#####################################
########## Set Parameters
#####################################
bh_dm_dir <- "/home2/sjseo/Bohun/SNParray/data/BH_DM/"
bh_nm_dir <- "/data/bohun_methylation/5566/메틸레이션_질본_정상_추가/KCHIP_836/"

dr_dir <- "/home2/sjseo/Bohun/SNParray/data/"

ncores <- 4

#####################################
########## Data Preparation
#####################################

bim_dm <- read.table(paste0(bh_dm_dir, "BH_DM_merged.bim"), stringsAsFactors = F)
bim_nm <- read.table(paste0(bh_nm_dir, "BH_836.bim"), stringsAsFactors = F)


##### Match SNP ID
upd <- mclapply(1:22, function(i) {
  #i <- 1
  chr_dm <- bim_dm %>% filter(V1 == i)
  chr_nm <- bim_nm %>% filter(V1 == i)
  
  chr_dm_upd <- chr_dm %>% left_join(., chr_nm, by = "V4") %>%
    mutate(V2.xx = ifelse((V5.x == V5.y & V6.x == V6.y), V2.y, V2.x))
  
  common_snp <- na.omit(chr_dm_upd$V2.xx)
  
  chr_dm_upd <- chr_dm_upd %>% 
    mutate(V2.xx = ifelse(is.na(V2.xx), V2.x, V2.xx)) %>%
    select(V1.x, V2.xx, V3.x, V4, V5.x, V6.x)
  
  return(list(chr_dm_upd, common_snp))
}, mc.cores = ncores)

bim_dm_upd <- do.call(rbind, lapply(upd, function(x) x[[1]]))
write.table(bim_dm_upd, paste0(bh_dm_dir, "BH_DM_merged.bim"), col.names = F, row.names = F, quote = F)

common_snps <- unlist(lapply(upd, function(x) x[[2]]))
write.table(common_snps, paste0(dr_dir, "01_02_Common_SNPs.txt"), sep = "\n", col.names = F, row.names = F, quote = F)


##### Merge DM and Normal data and select only common SNPs
merg_cmd <- paste("plink --bfile", paste0(bh_dm_dir, "BH_DM_merged"),
                  "--bmerge", paste0(bh_nm_dir, "BH_836"),
                  "--extract", paste0(dr_dir, "Common_SNPs.txt"),
                  "--geno 0.1 --hwe 0.001 --maf 0.05 --mind 0.05",
                  "--make-bed --out", paste0(dr_dir, "BH_merged"))
system(merge_cmd)

fam <- read.table(paste0(dr_dir,  "BH_merged.fam"), header = F, stringsAsFactors = F)
fam_upd <- fam %>% mutate(newID = ifelse(grepl("^GS", V1), V1, paste0(V1, "_", V2))) %>%
  mutate(V1 = newID, V2 = newID) %>% select(1:6)
write.table(fam_upd, paste0(dr_dir,  "BH_merged.fam"), col.names = F, row.names = F, quote = F)

##### Pruning
prune_cmd <- paste("plink --bfile", paste0(dr_dir, "BH_merged"),
                   "--indep-pairwise 500 5 0.5")
system(prune_cmd)

prune_out_cmd <- paste("plink --bfile", paste0(dr_dir, "BH_merged"),
                       "--extract", paste0(dr_dir, "plink.prune.in"),
                       "--make-bed --out", paste0(dr_dir, "BH_merged_pruned"))