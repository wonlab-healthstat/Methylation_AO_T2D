####################################
### Study : Bohun_Methylation
### Purpose of program : ID mathching (KoGES from bohun - KoGES from WonLab)
### Pre-required program : NA
### Programmed by : Sujin Seo
### Output : PLINK format file
### Draft Date : 2022-04-11
### Revision : NA
### Software (version) : R (v4.1.1)
#####################################

#####################################
########## Load Required Libraries
#####################################
library(data.table)
library(dplyr)
library(parallel)
library(readxl)

#####################################
########## Set Parameters
#####################################
raw_koges <- "/data2/KoGES_2021/KoGES_KARE_NC_CT_Kchip_Genotype/1.Genotyped/1.Raw/201125_12_Kchip_72297"
raw_bh_dm <- "/home2/sjseo/Bohun/SNParray/data/BH_DM/BH_DM_merged"
raw_bh_nm <- "/data/bohun_methylation/5566/메틸레이션_질본_정상_추가/KCHIP_836/BH_836"
raw_orange <- "/data/bohun_methylation/4171/multiomics/멀티오믹스를 이용한 근감소증 바이오마커발굴/6. 고엽제 당뇨 비고엽제 당뇨 한국인칩을 이용한 genetic risk score, polygenic risk score 비교/K-chip data/defoliant"

impute_koges <- "/data4/ardo/Data2/20201102_ImpuatationTest/NC_CT_KARE_Bohun_Merge/Merge_kareAffy_ncAffy_ctAffy_bhKchip_kncKchip/4.snpQC_genomafhwe/"
dr_dir <- "/home2/sjseo/Bohun/SNParray/derived_data/02_01_IDmatch/"

ncores <- 20

#####################################
########## Read files
#####################################

########## .bim files
koges <- fread(paste0(raw_koges, ".bim"), header = F, 
               stringsAsFactors = F, data.table = F)  # 467,088
bh_dm <- fread(paste0(raw_bh_dm, ".bim"), header = F, 
                stringsAsFactors = F, data.table = F) # 8,056,211
bh_nm <- fread(paste0(raw_bh_nm, ".bim"), header = F, 
               stringsAsFactors = F, data.table = F)  # 549,435

########## .fam files
orange_fam <- fread(paste0(raw_orange, ".fam"), header = F, 
                stringsAsFactors = F, data.table = F) # 804,105
impute_fam <- fread(paste0(impute_koges, "cleaned_chr1_Merge_kareAffy_ncAffy_ctAffy_bhKchip_kncKchip.fam"), header = F, 
                    stringsAsFactors = F, data.table = F)

########## other files
orange_match <- read_xlsx("/data/bohun_methylation/EPC array total information_rex.xlsx", sheet = 1)


########## Number of cores for parellel computing
ncores <- 20

#####################################
########## Extract common SNPs
#####################################

############# BH_DM + KoGES by chromosome/bp
##### Common SNPs
koges <- koges %>% mutate(chr_bp = paste0(V1, ":", V4))
bh_dm <- bh_dm %>% mutate(chr_bp = paste0(V1, ":", V4))

comb_dm <- koges %>% full_join(., bh_dm, by = "chr_bp") %>%
  mutate(Common = ifelse(V1.x == V1.y & V4.x == V4.y & V5.x == V5.y & V6.x == V6.y, "Common", "NotCommon"))

common_snp_dm <- comb_dm %>% filter(Common == "Common") #460,550
  
common_dm_koges <- common_snp_dm %>% select(V2.x)
common_dm_bh <- common_snp_dm %>% select(V2.y, V2.x)

write.table(common_dm_koges, paste0(dr_dir, "02_01_common_snp_dm_koges.txt"), col.names = F, row.names = F, quote = F)
write.table(common_dm_bh, paste0(dr_dir, "02_01_common_snp_dm_bh.txt"), col.names = F, row.names = F, quote = F)


##### Extract common SNPs
ext_common1 <- paste("plink --bfile", raw_koges, 
                    "--extract", paste0(dr_dir, "common_snp_dm_koges.txt"), 
                    "--make-bed --out", paste0(dr_dir, "KoGES_DM_common"))
system(ext_common1)

ext_common2 <- paste("plink --bfile", raw_bh_dm, 
                     "--extract", paste0(dr_dir, "common_snp_dm_bh.txt"), 
                     "--update-name", paste0(dr_dir, "common_snp_dm_bh.txt"), 
                     "--make-bed --out", paste0(dr_dir, "Bohun_DM_common"))
system(ext_common2)


##### Merge data
mrg_dm <- paste("plink --bfile", paste0(dr_dir, "KoGES_DM_common"),
                "--bmerge", paste0(dr_dir, "Bohun_DM_common"), 
                "--make-bed --out", paste0(dr_dir, "Merge_DM"))
system(mrg_dm)

##### Prunning
prune_dm <- paste("plink --bfile", paste0(dr_dir, "Merge_DM"), 
                  "--indep-pairwise 50 5 0.5", 
                  "--out", paste0(dr_dir, "Prune_DM"))
system(prune_dm)

prune_dm <- paste("plink --bfile", paste0(dr_dir, "Merge_DM"), 
                  "--extract", paste0(dr_dir, "Prune_DM.prune.in"), 
                  "--make-bed --out", paste0(dr_dir, "Pruned_DM"))
system(prune_dm)

############# BH_NM + KoGES by chromosome/bp
##### Common SNPs
koges <- koges %>% mutate(chr_bp = paste0(V1, ":", V4))
bh_nm <- bh_nm %>% mutate(chr_bp = paste0(V1, ":", V4))

comb_nm <- koges %>% full_join(., bh_nm, by = "chr_bp") %>%
  mutate(Common = ifelse(V1.x == V1.y & V4.x == V4.y & V5.x == V5.y & V6.x == V6.y, "Common", "NotCommon"))

common_snp_nm <- comb_nm %>% filter(Common == "Common") #454,311

common_nm_koges <- common_snp_nm %>% select(V2.x)
common_nm_bh <- common_snp_nm %>% select(V2.y, V2.x)

write.table(common_nm_koges, paste0(dr_dir, "common_snp_nm_koges.txt"), col.names = F, row.names = F, quote = F)
write.table(common_nm_bh, paste0(dr_dir, "common_snp_nm_bh.txt"), col.names = F, row.names = F, quote = F)


##### Extract common SNPs
ext_common1 <- paste("plink --bfile", raw_koges, 
                     "--extract", paste0(dr_dir, "common_snp_nm_koges.txt"), 
                     "--make-bed --out", paste0(dr_dir, "KoGES_nm_common"))
system(ext_common1)

ext_common2 <- paste("plink --bfile", raw_bh_nm, 
                     "--extract", paste0(dr_dir, "common_snp_nm_bh.txt"), 
                     "--update-name", paste0(dr_dir, "common_snp_nm_bh.txt"), 
                     "--make-bed --out", paste0(dr_dir, "Bohun_nm_common"))
system(ext_common2)


##### Merge data
mrg_nm <- paste("plink --bfile", paste0(dr_dir, "KoGES_nm_common"),
                "--bmerge", paste0(dr_dir, "Bohun_nm_common"), 
                "--make-bed --out", paste0(dr_dir, "Merge_nm"))
system(mrg_nm)

##### Prunning
prune_nm <- paste("plink --bfile", paste0(dr_dir, "Merge_nm"), 
                  "--indep-pairwise 50 5 0.5", 
                  "--out", paste0(dr_dir, "Prune_nm"))
system(prune_nm)


prune_nm <- paste("plink --bfile", paste0(dr_dir, "Merge_nm"), 
                  "--extract", paste0(dr_dir, "Prune_nm.prune.in"), 
                  "--make-bed --out", paste0(dr_dir, "Pruned_nm"))
system(prune_nm)




#####################################
########## IBS
#####################################
##### Calculate IBS
ibs_dm <- paste("plink --bfile", paste0(dr_dir, "Pruned_DM"), 
                "--cluster --matrix --genome --min 0.5", 
                "--out", paste0(dr_dir, "IBS_dm"))
system(ibs_dm)

ibs_nm <- paste("plink --bfile", paste0(dr_dir, "Pruned_nm"), 
                "--cluster --matrix --genome --min 0.5", 
                "--out", paste0(dr_dir, "IBS_nm"))
system(ibs_nm)



##### Set matching threshold
ibs_dm <- fread(file.path(dr_dir, "IBS_dm.genome"), head=T, stringsAsFactors = F, data.table=F) 
ibs_nm <- fread(file.path(dr_dir, "IBS_nm.genome"), head=T, stringsAsFactors = F, data.table=F) # all are PI_HAT == 1

### Decide threshold based on graph
summary(ibs_dm$PI_HAT)
thr <- seq(98, 100, by = 0.01)
thr_n <- sapply(thr, function(i) {
  n <- sum(ibs_dm$PI_HAT >= i/100)
}) 
plot(thr, thr_n/nrow(ibs_dm), xlab = "Threshold", ylab = "Matching rate", ylim = c(0, 1))
unique(thr_n/nrow(ibs_dm))
abline(h = 0.98, col = "red")
thr_dm <- 0.98


#####################################
########## ID list
#####################################

id_list <- rbind(
  ##### DM with AO
  orange_match %>% filter(ID %in% impute_fam$V2) %>%
    mutate(Source = "Bohun", Pheno = "DMwithAO") %>% select(Source, Pheno, Methyl = Sample_ID, Geno_Bohun = ID, Geno_Wonlab = ID),
  ##### DM without AO
  ibs_dm %>% mutate(Source = "KCDC", Pheno = "DMwithoutAO") %>% select(Source, Pheno, Methyl = IID1, Geno_Bohun = IID1, Geno_Wonlab = IID2),
  ##### Normal
  ibs_nm %>% mutate(Source = "KCDC", Pheno = "Normal") %>% select(Source, Pheno, Methyl = IID1, Geno_Bohun = IID1, Geno_Wonlab = IID2)
)
write.table(id_list, 
            paste0(dr_dir, "Matched_ID.txt"), sep = "\t", col.names = T, row.names = F, quote = F)


##### Idential with the above but for PLINK 
id_list_plink <- rbind(ibs_dm %>% select(V1 = FID2, V2 = IID2),                         # DM
                 ibs_nm %>% select(V1 = FID2, V2 = IID2),                               # Normal 
                 impute_fam %>% filter(V2 %in% orange_match$ID) %>% select(V1, V2))     # DM wi AO 
write.table(id_list_plink, 
            paste0(dr_dir, "Matched_ID_plink.txt"), sep = "\t", col.names = F, row.names = F, quote = F)


#####################################
########## Extract matched subject from 70K data
#####################################
id_list <- read.table(paste0(dr_dir, "Matched_ID.txt"), sep = "\t", header = T, stringsAsFactors = F)

mclapply(1:22, function(i) {
  ext_subj_dm <- paste("plink --bfile", paste0(impute_koges, "cleaned_chr", i, "_Merge_kareAffy_ncAffy_ctAffy_bhKchip_kncKchip"), 
                       "--keep", paste0(dr_dir, "Matched_ID_plink.txt"), 
                       "--make-bed --out", paste0(dr_dir, "Matched_chr", i))
  system(ext_subj_dm)
  
  org_fam <- read.table(paste0(dr_dir, "Matched_chr", i, ".fam")) %>% 
    left_join(., id_list, by = c("V2"="Geno_Wonlab")) %>%
    select(V1 = Methyl, V2 = Methyl, V3, V4, V5, V6)
  
  write.table(org_fam, paste0(dr_dir, "Matched_chr", i, ".fam"),
              col.names = F, row.names = F, quote = F)
}, mc.cores = ncores)


#####################################
########## Merged extracted files
#####################################
file_list <- unlist(lapply(2:22, function(i) paste0(paste0(dr_dir, "Matched_chr", i), 
                                                   c(".bed", ".bim", ".fam"), 
                                                   collapse = " ")))

file_list <- paste0(dr_dir, "Matched_chr", 2:22)
write.table(file_list, paste0(dr_dir, "File_list.txt"), sep = "\n", col.names = F, row.names = F, quote = F)

mrg <- paste("plink --bfile", paste0(dr_dir, "Matched_chr1"),
             "--merge-list", paste0(dr_dir, "File_list.txt"), 
             "--make-bed", 
             "--out", paste0(dr_dir, "../02_01_IDmatched"))
system(mrg)
