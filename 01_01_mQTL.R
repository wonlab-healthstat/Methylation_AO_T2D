#####################################
### Study : Bohun-Methylation
### Purpose of prgram : Methylation-SNP QTL
### Pre-required program : NA
### Output files : 
### Programmed by : Sujin Seo
### Draft Date : 2021-02-25
### Revision : NA
### Software (version) : R (v4.1.1)
#####################################

#####################################
##### Load Required Libraries
#####################################

library(MatrixEQTL)
library(dplyr)
library(data.table)
library(parallel)


#####################################
##### Set Parameters
#####################################

raw_snp_dir <- "/home2/sjseo/Bohun/SNParray/derived_data/"
raw_meth_dir <- "/home2/sjseo/Bohun/Methylation/derived_data/"
raw_meth_dir2 <- "/data/bohun_methylation/4171/multiomics/멀티오믹스를 이용한 근감소증 바이오마커발굴/5. 고엽제 당뇨 EWAS 비고엽제당뇨 EWAS Epic 850K/2차 당뇨 메틸레이션 72명"
raw_meth_dir3 <- "/data/bohun_methylation/4171/multiomics/멀티오믹스를 이용한 근감소증 바이오마커발굴/5. 고엽제 당뇨 EWAS 비고엽제당뇨 EWAS Epic 850K/1차 당뇨 메틸레이션 53명/EPIC array 1차 53명"
dr_dir <- "/home2/sjseo/Bohun/mQTL/derived_data/"

ncores <- 20


#################################################
########## Data preparation
#################################################

########## Methylation 
##### Match sample ID with barcode ID
meth_smp1 <- read.table(file.path(raw_meth_dir2, "Sample.Table.txt"), sep = "\t", header = T, stringsAsFactors = F)
meth_smp2 <- read.table(file.path(raw_meth_dir3, "Sample.Table 1.txt"), sep = "\t", header = T, stringsAsFactors = F)
meth_smp3 <- read.table(file.path(raw_meth_dir3, "Sample.Table 2.txt"), sep = "\t", header = T, stringsAsFactors = F)
meth_smp <- rbind(meth_smp1, meth_smp2, meth_smp3) %>%
  mutate(Barcode = paste(Sentrix.Barcode, Sample.Section, sep = "_")) %>% 
  select(Barcode, Sample.ID)


##### Methylation data
manifest <- fread(paste0(raw_meth_dir, "manifest.csv"))
load(paste0(raw_meth_dir, "03_02_QCBatch_M_combat.RData"))

barcode_id <- data.frame(Methyl_ID = colnames(batch_cntrl1)) %>%
  left_join(meth_smp, by = c("Methyl_ID" = "Barcode")) %>%
  mutate(Sample.ID = ifelse(is.na(Sample.ID), Methyl_ID, Sample.ID))
colnames(batch_cntrl1) <- barcode_id$Sample.ID 
save(batch_cntrl1, file = paste0(dr_dir, "01_01_mQTL_methyl_dat.RData"))


### Target CpG
sig_cpg <- read.table(paste0(raw_meth_dir, "04_01_Limma_Subgroup_Table_Nested.csv"), sep = ",", header = T, stringsAsFactors = F)  
a <- sig_cpg %>% filter(Type == "defoliant_DM", significance != "Not Significant") %>%
  left_join(., manifest %>% select(probe_id, chr, mapinfo), by = c("CpG" = "probe_id"))
b <- sig_cpg %>% filter(Type == "DM", significance != "Not Significant") %>%
  left_join(., manifest %>% select(probe_id, chr, mapinfo), by = c("CpG" = "probe_id"))
target <- intersect(a$CpG, b$CpG)
target_chr <- na.omit(as.numeric(unique(b$chr)))
target_meth <- batch_cntrl1[b$CpG, ]
rm("batch_cntrl1")


########## SNP data
fam <- read.table(file.path(raw_snp_dir,  "02_01_IDmatched.fam"), header = F, stringsAsFactors = F)

##### Extract common subjects (QCed methyl vs Matched genotype)
pheno <- read.table(paste0(raw_meth_dir, "QC_Pheno.csv"), sep = ",",
                    header=TRUE, stringsAsFactors=FALSE) %>%
  filter(QC_final == "PASS", PC_out_combat == "Include") %>%
  left_join(., meth_smp, by = c("Sample_ID" = "Barcode")) %>%
  mutate(Sample_ID = ifelse(is.na(Sample.ID), Sample_ID, Sample.ID)) %>%
  select(-Sample.ID) #1,249

write.table(pheno, paste0(raw_meth_dir, "QC_Pheno.csv"), sep = ",", 
            col.name = T, row.names = F, quote = F)

both_iid <- intersect(fam$V2, pheno$Sample_ID) #1,106

write.table(fam %>% filter(V2 %in% both_iid) %>% select(V1, V2), 
            paste0(dr_dir, "01_01_Methyl_Common_ID.txt"), sep = "\t", col.names = F, row.names = F, quote = F)

##### Pruning
cmd_prune <- paste("plink --bfile", paste0(raw_snp_dir, "02_01_IDmatched"), 
                   "--keep", paste0(dr_dir, "01_01_Methyl_Common_ID.txt"),
                   "--indep-pairwise 500 5 0.5", 
                   "--out", paste0(raw_snp_dir, "Prune"))
system(cmd_prune)
# pruned data were deleted after running cod 01_01_mQTL.R

cmd_ext_prune <- paste("plink --bfile", paste0(raw_snp_dir, "IDmatched"),
                       "--keep", paste0(dr_dir, "Methyl_Common_ID.txt"),
                       "--extract", paste0(raw_snp_dir, "Prune.prune.in"),
                       "--make-bed", 
                       "--out", paste0(raw_snp_dir, "IDmatched_Prune"))
system(cmd_ext_prune)
# 2,144,245 SNPs -> 121,645 SNPs for 1,106 individuals
# pruned data were deleted after running cod 01_01_mQTL.R

### SNP information
snppose <- read.table(paste0(raw_snp_dir, "IDmatched_Prune.bim"), sep = "\t", header = F, stringsAsFactors = F)
# pruned data were deleted after running cod 01_01_mQTL.R

########## Covariate file for eQTL
smp_idx <- sapply(both_iid, function(x) which(pheno$Sample_ID == x))
pheno_both <- pheno[smp_idx,] %>% select(Sample_ID, SEX, DM) %>%
  mutate(SEX = as.numeric(as.factor(SEX)), DM = as.numeric(DM)) %>%
  tibble::remove_rownames(.) %>%
  tibble::column_to_rownames(var = "Sample_ID") %>% t()
write.table(pheno_both, paste0(dr_dir, "01_01_Pheno_QTL.txt"), sep = "\t", col.names = T, row.names = T, quote = F)



#################################################
########## eQTL
#################################################

nTest <- data.frame()
for(i in target_chr) {
  recode_cmd <- paste("plink --bfile", paste0(raw_snp_dir, "IDmatched_Prune"), 
                      "--chr", i, 
                      "--recode A-transpose --out", paste0(dr_dir, "01_01_Prune/Prune_chr", i))
  system(recode_cmd)
  
  snp_dat <- fread(paste0(dr_dir, "01_01_Prune/Prune_chr", i, ".traw"))
  snp_dat_upd <- snp_dat %>% tibble::column_to_rownames(var = "SNP") %>% select(-c(1:5))
  snp_iid <- colnames(snp_dat_upd)
  snp_iid2 <- sapply(snp_iid, function(x) {
    idx <- (nchar(x)-1)/2
    substr(x, 1, idx)
  })
  colnames(snp_dat_upd) <- snp_iid2

  snppos_chr <- snp_dat_upd %>% tibble::rownames_to_column("snpid") %>% select(snpid) %>%
    left_join(., snppose, by = c("snpid" = "V2")) %>% 
    select(snpid, chr = V1, pos = V4)
  nSNP <- nrow(snppos_chr)
  
  write.table(snp_dat_upd, paste0(dr_dir, "01_01_SNP/SNP_chr", i, ".txt"), sep = "\t", col.names = T, row.names = T, quote = F)
  system(paste("rm", paste0(dr_dir, "01_01_Prune/Prune_chr", i, ".traw")))
  rm(list = c("snp_dat", "snp_dat_upd"))
  
  ##### Methylation
  chr_prob <- b[b$chr == i, "CpG"]
  id_idx <- unlist(sapply(snp_iid2, function(x) which(colnames(target_meth) == x)))
  arrange_expr <- target_meth[rownames(target_meth) %in% chr_prob, id_idx]
  genepos_chr <- arrange_expr %>% tibble::rownames_to_column("geneid") %>% select(geneid) %>%
    left_join(., manifest, by = c("geneid" = "probe_id")) %>% 
    select(geneid, chr, left = mapinfo, right = mapinfo)
  nCpG <- nrow(genepos_chr)
  write.table(arrange_expr, paste0(dr_dir, "01_01_Methylation/Methylation_chr", i, ".txt"), sep = "\t", col.names = T, row.names = T, quote = F)
  rm("arrange_expr")
  
  
  #################################################
  ########## eQTL
  #################################################
  
  ##### Set parameters
  pvOutputThreshold = 1
  pvOutputThreshold.cis = 1
  errorCovariance = numeric()
  useModel = modelLINEAR
  
  ##### SNP
  SNP_file_name = paste0(dr_dir, "01_01_SNP/SNP_chr", i, ".txt")
  snps = SlicedData$new()
  snps$fileDelimiter = "\t"     # the TAB character
  snps$fileOmitCharacters = "NA" # denote missing values;
  snps$fileSkipRows = 1          # one row of column labels
  snps$fileSkipColumns = 1       # six columns of row labels
  snps$fileSliceSize = 2000     # read file in pieces of 2,000 rows
  snps$LoadFile(SNP_file_name)
  
  
  ##### Methylation
  expression_file_name = paste0(dr_dir, "01_01_Methylation/Methylation_chr", i, ".txt")
  gene = SlicedData$new()
  gene$fileDelimiter = "\t"      # the TAB character
  gene$fileOmitCharacters = "NA" # denote missing values;
  gene$fileSkipRows = 1          # one row of column labels
  gene$fileSkipColumns = 1      # one column of row labels
  gene$fileSliceSize = 2000      # read file in pieces of 2,000 rows
  gene$LoadFile(expression_file_name)
  
  ##### Covariates
  covariates_file_name = paste0(dr_dir, "01_01_Pheno_QTL.txt")
  cvrt = SlicedData$new()
  cvrt$fileDelimiter = "\t"      # the TAB character
  cvrt$fileOmitCharacters = "NA" # denote missing values;
  cvrt$fileSkipRows = 1          # one row of column labels
  cvrt$fileSkipColumns = 1      # one column of row labels
  cvrt$fileSliceSize = 2000      # read file in pieces of 2,000 rows
  cvrt$LoadFile(covariates_file_name)
  
  ##### Run MatrixEQTL
  out_file_cis = paste0(dr_dir, "01_01_mQTL/mQTL_chr", i, "_cis.txt")
  out_file_trans = paste0(dr_dir, "01_01_mQTL/mQTL_chr", i, "_trans.txt")
  
  me = Matrix_eQTL_main(
    snps = snps, snpspos = snppos_chr, 
    gene = gene, genepos = genepos_chr,
    cvrt = cvrt,
    output_file_name = out_file_trans,
    output_file_name.cis = out_file_cis,
    pvOutputThreshold = pvOutputThreshold,
    pvOutputThreshold.cis = pvOutputThreshold.cis,
    cisDist = 500000,
    useModel = useModel, 
    errorCovariance = errorCovariance, 
    verbose = TRUE,
    pvalue.hist = "qqplot", 
    min.pv.by.genesnp = FALSE,
    noFDRsaveMemory = FALSE);
  
  png(paste0(dr_dir, "01_01_mQTL_QQplot/mQTL_QQplot_chr", i, ".png"))
  plot(me)
  dev.off()
  
  nTest <- rbind(nTest, 
                 data.frame(Chr = i, nSNP, nCpG, nTests = me$all$ntests, nTests_trans = me$trans$ntests, nTests_cis = me$cis$ntests))
  cat(paste('Chr', i, 'Detected local eQTLs:', me$cis$ntests, '\n'));
  cat(paste('Chr', i, 'Detected distant eQTLs:', me$trans$ntests, '\n'));
  
  rm("me")
}
write.table(nTest %>% arrange(Chr), paste0(dr_dir, "01_01_nTest.txt"), sep = "\t", col.names = T, row.names = F, quote = F)


