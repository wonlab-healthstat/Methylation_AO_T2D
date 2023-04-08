#####################################
### Study : Bohun-Methylation
### Purpose of prgram : Logistic regression : DM~SNP for the significant SNPS in mQTL
### Pre-required program : NA
### Output files : 
### Programmed by : Sujin Seo
### Draft Date : 2021-03-07
### Revision : NA
### Software (version) : R (v4.1.1)
#####################################


#####################################
##### Load Required Libraries
#####################################

library(dplyr)
library(data.table)
library(ggplot2)
library(parallel)
#library(systemfit)
#library(ivtools)
library(ivreg)
library(tibble)
library("cowplot")
source("/home2/sjseo/Tools/R/library/qqunif.plot.R")


#####################################
##### Set Parameters
#####################################

dr_dir <- "/home2/sjseo/Bohun/mQTL/derived_data/"
raw_snp_dir <- "/home2/sjseo/Bohun/SNParray/derived_data/"
raw_meth_dir <- "/home2/sjseo/Bohun/Methylation/derived_data/"
raw_meth_dir2 <- "/data/bohun_methylation/4171/multiomics/멀티오믹스를 이용한 근감소증 바이오마커발굴/5. 고엽제 당뇨 EWAS 비고엽제당뇨 EWAS Epic 850K/2차 당뇨 메틸레이션 72명"
raw_meth_dir3 <- "/data/bohun_methylation/4171/multiomics/멀티오믹스를 이용한 근감소증 바이오마커발굴/5. 고엽제 당뇨 EWAS 비고엽제당뇨 EWAS Epic 850K/1차 당뇨 메틸레이션 53명/EPIC array 1차 53명"

ncores <- detectCores()



#########################################################
##### Mendelian Randomization (Single)
#########################################################


### Methylation data
#load(paste0(raw_meth_dir, "QCBatch_M_combat.RData"))
manifest <- fread(paste0(raw_meth_dir, "manifest.csv"))



##### Match sample ID with barcode ID
meth_smp1 <- read.table(file.path(raw_meth_dir2, "Sample.Table.txt"), sep = "\t", header = T, stringsAsFactors = F)
meth_smp2 <- read.table(file.path(raw_meth_dir3, "Sample.Table 1.txt"), sep = "\t", header = T, stringsAsFactors = F)
meth_smp3 <- read.table(file.path(raw_meth_dir3, "Sample.Table 2.txt"), sep = "\t", header = T, stringsAsFactors = F)
meth_smp <- rbind(meth_smp1, meth_smp2, meth_smp3) %>%
  mutate(Barcode = paste(Sentrix.Barcode, Sample.Section, sep = "_")) %>% 
  dplyr::select(Barcode, Sample.ID)

pheno <- read.table(paste0(raw_meth_dir, "04_01_Phenotypes_1249.csv"), sep = ",",
                    header=TRUE, stringsAsFactors=FALSE) %>%
  #filter(QC_final == "PASS", PC_out_combat == "Include") %>%
  left_join(., meth_smp, by = c("Sample_ID" = "Barcode")) %>%
  mutate(Sample_ID = ifelse(is.na(Sample.ID), Sample_ID, Sample.ID), 
         defoliant = ifelse(CASE2 == "Case", 1, 0), 
         CASE = ifelse(defoliant == 1, "Case", ifelse(DM == 1, "Control", "Normal"))) %>%
  dplyr::select(-Sample.ID)

### Target CpG
load(paste0(dr_dir, "mQTL_methyl_dat.RData"))
sig_cpg <- read.table(paste0(raw_meth_dir, "04_01_Limma_Subgroup_Table_Nested.csv"), sep = ",", header = T, stringsAsFactors = F)  
a <- sig_cpg %>% filter(Type == "defoliant_DM", significance != "Not Significant") %>%
  left_join(., manifest %>% dplyr::select(probe_id, chr, mapinfo), by = c("CpG" = "probe_id"))
b <- sig_cpg %>% filter(Type == "DM", significance != "Not Significant") %>%
  left_join(., manifest %>% dplyr::select(probe_id, chr, mapinfo), by = c("CpG" = "probe_id"))
target <- intersect(a$CpG, b$CpG)
target_chr <- na.omit(as.numeric(unique(b$chr)))
target_meth <- batch_cntrl1[b$CpG, ]
#meth <- batch_cntrl1 ; rm(batch_cntrl1)
rm("batch_cntrl1")


  
#### GWAS results for reverse causation
#prs <- read.table("/home2/sjseo/Bohun/SNParray/derived_data/03_03_PRS_Wonlab/IDmatched.LDpred.inf.prs", header = T, stringsAsFactors = F)

gwas_res <- read.table("/home2/sjseo/Bohun/SNParray/derived_data/03_03_PRS_Wonlab/DM_discovery.sumstats",
                       sep = "\t", stringsAsFactors = F, header= T) %>%
  mutate(adj_pval = p.adjust(gwas_res$P, method = "BH"))

res_cis_rs <- read.table(paste0(dr_dir, "01_02_mQTL_cis_rsID.txt"), header = T,stringsAsFactors = F)

snp_dat <- lapply(1:22, function(i) {
  snp_dat <- read.table(paste0(dr_dir, "01_01_SNP/SNP_chr", i, ".txt"))
  return(snp_dat)
})
snp_dat <- do.call(rbind, snp_dat)


snp_list <- gwas_snps <- gwas_res %>% 
  filter(adj_pval < 0.05, 
         !(SNP %in% res_cis_rs[res_cis_rs$padj < 0.05, "SNP"]),
         SNP %in% rownames(snp_dat)) 


#final_target <- snp_list %>% filter(gene %in% target, padj < 0.05)
#sem <- lapply(1:nrow(snp_list), function(i) {
sem <- mclapply(1:nrow(b), function(i) {

  target_cpg <- b[i,"CpG"]
  tmp_chr <- b[i, 'chr']
  tmp_mapinfo <- b[i, 'mapinfo']
  
  ### Reverse Causation
  # IV : PRS
  # Y (Methylation)
  tmp_meth_dat <- t(target_meth[target_cpg, ]) %>% data.frame() %>%
    rownames_to_column("IID")
  # X (T2D) + Combine above
  tmp_dat <- prs %>% left_join(., tmp_meth_dat, by = "IID") %>%
    left_join(., pheno %>% dplyr::select(IID = Sample_ID, SEX, AGE, CASE2), by = "IID")
  tmp_dat <- tmp_dat[complete.cases(tmp_dat),] 
  colnames(tmp_dat)[3:4] <- c("PRS", "Meth")
  tmp_dat$CASE2 <- factor(tmp_dat$CASE2, levels = c("Control", "Case"))
  
  if(all(table(tmp_dat$CASE2)) > 0) {
    ### Reverse causation
    ive_reg <- ivreg(Meth ~  CASE2 + SEX + AGE | PRS + SEX + AGE, data = tmp_dat)
    ive_coef <- (coef(summary(ive_reg)) %>% data.frame())[2,]
    
    ive_coef_res2 <- data.frame(CHR = tmp_chr, Meth = target_cpg, MAPINFO = tmp_mapinfo,
                                Pheno_Meth_Est = ive_coef[, 1], Pheno_Meth_SE = ive_coef[, 2], 
                                Pheno_Meth_Stat = ive_coef[, 3], Pheno_Meth_Pval = ive_coef[, 4])
    
    return(ive_coef_res2)
  }
}, mc.cores = ncores)

#})

# https://cran.r-project.org/web/packages/ivreg/vignettes/ivreg.html
# https://github.com/alicerosecarter/MediationMR/blob/master/example_R.R

sem_reverse <- do.call(rbind, lapply(sem, function(x) x)) %>%
  mutate(Pheno_Meth_Adj_Pval = p.adjust(Pheno_Meth_Pval)) %>% 
  arrange(Pheno_Meth_Pval)

sem_reverse %>% filter(Pheno_Meth_Adj_Pval < 0.05)


write.table(sem_reverse, paste0(dr_dir, "02_02_MR_reverse_single_0927.txt"), sep = "\t", 
            col.names = T, row.names = F, quote = F)


