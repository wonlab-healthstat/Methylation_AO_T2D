#####################################
### Study : Bohun-Methylation
### Purpose of prgram : Quality Control Methylation Data
### Pre-required program : NA
### Output files : Start with "QC"
### Programmed by : Sujin Seo
### Draft Date : 2021-10-07
### Revision : NA
### Software (version) : R (v3.6.3)
#####################################

#####################################
##### Load Required Libraries
#####################################
library(ewastools)
library(meffil)
library(dplyr)
library(readr)
library(ggplot2)
library(svd)
library(parallel)
library(ChAMP)
library(wateRmelon)


#####################################
##### Set Parameters
#####################################

##### Directories
dr_dir <- "/home2/sjseo/Bohun/Methylation/derived_data/"
meth_dir <- paste0("/home2/sjseo/Bohun/Methylation/derived_data/", "meth_rdata/")

##### Number of cores
ncores <- 20

##### Output file
out_file <- "03_01_QC"


#####################################
##### Defined Function Code for QC
#####################################

QC_report <- function(meth) {
  
  manifest <- meth$manifest
  
  #####################################
  ##### Preprocess
  #####################################
  
  ##### Pre1. Detection p-values (Distinguish signals from noise)
  ##### Pre2. Dye-bias correction
  ##### Pre3. Do not normalize to prevent removal of genuine biological signals
  p_thr <- 0.01
  beta <- ewastools::detectionP(meth) %>% 
    ewastools::mask(p_thr) %>% 
    correct_dye_bias %>%
    dont_normalize
  #beta_df <- as.data.frame(beta) %>% mutate(probe = rownames(beta))
  # =>  The above code results in beta-values.
  
  
  #####################################
  ##### Quality Control
  #####################################
  
  ##### Probes
  ### QC_P1. Remove probes where # of beads are to small 
  #       (threshold : beadcount < 3 in 5% of samples)
  
  ##### Samples
  ### QC_S1. BeadArray Control Reports
  ctrls <- control_metrics(meth)
  
  ctrls_df <- lapply(1:17, function(i) {
    tmp_dat <- data.frame(ctrl_nm = names(ctrls)[i],
                          y = ctrls[[i]],
                          Sample_ID = meth$meta$sample_id,
                          threshold = attributes(ctrls[[i]])$threshold,
                          stringsAsFactors = F)
    return(tmp_dat)
  }) 
  ctrls_df <- do.call("rbind", ctrls_df) 
  
  
  ### QC_S2. Sex Mismatches
  xy_int <- check_sex(meth) # Normalized average intensity of X, Y chr
  pred_sex <- data.frame(Sample_ID = meth$meta$sample_id, 
                         Pred_Sex_X = xy_int$X,
                         Pred_Sex_Y = xy_int$Y,
                         stringsAsFactors = F)
  
  
  ### QC_S3. Genotype Calling
  snps <- meth$manifest[probe_type=="rs", "probe_id"] # Row index of SNP probes
  snps_index <- which(rownames(beta) == snps)
  
  genotypes <- call_genotypes(beta[snps_index, ], learn = FALSE)
  
  snp_log_odds <-  genotypes$outliers/(1 - genotypes$outliers)
  snp_outliers <- mean(log2(snp_log_odds), na.rm = TRUE)
  
  ### Return QC tables
  res <- list(Sample_ID = meth$meta$sample_id,  
              manifest = manifest, 
              beta = beta, 
              ctrl = ctrls_df,
              pred_sex = pred_sex, 
              snp_outlier = snp_outliers)
  return(res)
}



#####################################
##### Run QC
#####################################

files <- list.files(meth_dir, full.names = T)

#### Return QC
QC_tab <- mclapply(files[1:3], function(f) {
  nm <- strsplit(strsplit(f, paste0(dr_dir, "meth_rdata//"))[[1]][2], 
                 ".RData")[[1]][1]
  print(paste(f, "START"))
  load(f)
  QC_tab <- QC_report(tmp_dat)
  save(QC_tab, file = paste0(dr_dir, out_file, "_tab/", nm, ".RData"))
  rm(list = c("QC_tab", "tmp_dat"))
  return("SUCCESS")
}, mc.cores = ncores)
# Delete QC_tab folder after run 03_01_QC.R for memory saving

#####################################
##### QC Report
#####################################

##### Betas
QC_files <- list.files(paste0(dr_dir, out_file, "_tab"), full.names = T)

beta_abs <- mclapply(QC_files, function(f) {
  load(f)
  beta <- QC_tab$beta
  smp <- QC_tab$Sample_ID
  return(list(beta = beta, smp = smp))
}, mc.cores = ncores)

smp <- unlist(lapply(beta_abs, function(x) x[[2]]))
beta <- do.call("cbind",lapply(beta_abs, function(x) x[[1]]))
write.table(beta, paste0(dr_dir, out_file, "_beta_raw.csv"), col.names = T, row.names = T)

##### Pheno 
pheno <- read.table(paste0(dr_dir, "02_01_Read_idat_Pheno_Anal.csv"), sep = ",",
                    header = T, stringsAsFactors=FALSE) %>%
  filter(Sample_ID %in% smp) %>%
  mutate(CASE = factor(CASE, labels = c("Control", "Case")),
         CASE2 = ifelse(CASE == "Case", "Case",
                        ifelse(DM == 1, "Control", "Normal")),
         SEX = ifelse(SEX == 2, "F","M"))


##### QC_S1. BeadArray Control Reports
ctrls_df <- do.call("rbind", 
                    mclapply(QC_files, function(f) {
                      load(f)
                      ctrls <- QC_tab$ctrl
                      return(ctrls)
                    }, mc.cores = ncores)) %>%
  left_join(., pheno, by = "Sample_ID") %>% 
  mutate(QC = ifelse(y < threshold | is.na(y), "FAIL", "PASS"))

ctrls_fail <- lapply(smp, function(x) {
  tmp_dat <- ctrls_df %>% filter(Sample_ID == x)
  tmp_qc <- unique(tmp_dat$QC)
  
  if ("FAIL" %in% tmp_qc) {
    tmp_out <- data.frame(Sample_ID = x, QC = "FAIL")
  } else {
    tmp_out <- data.frame(Sample_ID = x, QC = "PASS")
  }
  return(tmp_out)
})
ctrls_fail <- do.call("rbind", ctrls_fail)
pheno <- pheno %>% left_join(., ctrls_fail %>% select(Sample_ID, QC_ctrl = QC))

##### QC_S2. Sex Mismatches
pred_sex <- do.call("rbind", 
                    mclapply(QC_files, function(f) {
                      load(f)
                      pred_sex <- QC_tab$pred_sex[, c("Sample_ID", "Pred_Sex_X", "Pred_Sex_Y")]
                      return(pred_sex)
                    }, mc.cores = ncores)) %>%
  left_join(., pheno, by = "Sample_ID") 
pred_sex <- pred_sex %>% 
  mutate(Pred_Sex = predict_sex(Pred_Sex_X, Pred_Sex_Y, 
                                which(SEX == "M"), which(SEX == "F")),
         QC = ifelse(toupper(Pred_Sex) == toupper(SEX), "PASS", "FAIL"))
pheno <- pheno %>% left_join(., pred_sex %>% select(Sample_ID, QC_sex = QC))


##### QC_S3. Genotype Calling and Outliers
snp_outlier <- unlist(mclapply(QC_files, function(f) {
  load(f)
  snp_out <- QC_tab$snp_outlier
  return(snp_out)
}, mc.cores = ncores))

snp_outlier <- data.frame(Sample_ID = smp, 
                          SNP_out = snp_outlier) %>%
  left_join(., pheno, by = "Sample_ID") %>%
  mutate(QC = ifelse(SNP_out <= -4, "PASS", "FAIL"))
pheno <- pheno %>% left_join(., snp_outlier %>% select(Sample_ID, QC_snpout = QC))


##### QC_S4. Estimate Cell Composition
LC <- meffil.estimate.cell.counts.from.betas(as.matrix(beta), 
                                             cell.type.reference = "andrews and bakulski cord blood",
                                             verbose=F)
write.csv(LC, paste0(dr_dir, out_file, "_LC.csv"), row.names = T, col.names = T)

lc_fract <- do.call("rbind", 
                    lapply(colnames(LC), function(x) {
                      data.frame(Sample_ID = rownames(LC), 
                                 Cell = x,
                                 fraction = LC[, x]*100)  
                    })) %>% left_join(., pheno, by = "Sample_ID")

### By manually review the graph, 
### Exclude a subject who has >30 on NK and >9 on nRBC
lc_fract_fail <- lc_fract %>% 
  filter((Cell == "NK" & fraction > 30) | (Cell == "nRBC" & fraction > 9))
pheno <- pheno %>% 
  mutate(QC_lc = ifelse(Sample_ID %in% lc_fract_fail$Sample_ID, "FAIL", "PASS"))


##### QC_S5. PCA
set.seed(2021)

# Exclude X&Y chromosome
load(QC_files[1])
manifest <- QC_tab$manifest ; rm("QC_tab")

chrXY <- manifest[manifest$chr %in% c("X","Y") & 
                    manifest$probe_type!="rs", "probe_id"]
chrXY_index <- rownames(beta) %in% chrXY
pcs <-  beta[!chrXY_index, ]

pcs <- t(na.omit(pcs - rowMeans(pcs)))
pc2 <- trlan.svd(pcs, neig=2) # compute the first two principal components
pc2_u <- data.frame(Sample_ID = rownames(pcs), 
                    PC1 = pc2$u[, 1],
                    PC2 = pc2$u[, 2],
                    stringsAsFactors = F) %>%
  left_join(., pheno, by = "Sample_ID")


##### Normalize Beta Values
data(probeInfoALL.epic.lv)
design.v <- as.numeric(lapply(probeInfoALL.lv, function(x) x)$Design[match(rownames(beta), 
                                                                           probeInfoALL.lv$probeID)])
if (min(beta, na.rm = TRUE) == 0) beta[beta == 0] <- 1e-06

norm_beta <- mclapply(1:ncol(beta), function(x) {
  if(x %% 10 == 0) print(paste(x, "in", ncol(beta), "START"))
  
  norm_beta <- BMIQ(beta[, x], design.v, 
                    sampleID = colnames(beta)[x], 
                    plots = FALSE, pri = FALSE)$nbeta
  
  if(x %% 10 == 0) print(paste(x, "in", ncol(beta), "END"))
  
  return(norm_beta)
}, mc.cores = ncores)
norm_beta <- do.call("cbind", norm_beta)
rownames(norm_beta) <- rownames(beta)
colnames(norm_beta) <- colnames(beta)       
write.table(norm_beta, paste0(dr_dir, out_file, "_beta_norm.csv"), col.names = T, row.names = T)


##### QC_P1. Filter probes with >3% missing
prob_miss <- unlist(mclapply(1:nrow(norm_beta), 
                             function(i) sum(is.na(norm_beta[i, ])),
                             mc.cores = ncores))
prob_out_idx <- which(prob_miss/ncol(norm_beta)*100 >= 3) 

##### QC_P2. Filter out the non-CG probes
prob_noncg <- manifest[manifest$probe_type != "cg", "probe_id"]
prob_noncg_idx <- which(rownames(norm_beta) %in% data.frame(prob_noncg)$probe_id)


##### Re-calculate PC scores with filtered data
pheno <- pheno %>% mutate(QC_final = ifelse((QC_ctrl == "FAIL" | QC_sex == "FAIL" | QC_snpout == "FAIL" | QC_lc == "FAIL"),
                                            "FAIL", "PASS"))
fail_smp <- pheno[pheno$QC_final == "FAIL", "Sample_ID"]

smp_out_idx <- which(colnames(norm_beta) %in% pheno[pheno$QC_final == "FAIL", "Sample_ID"])
beta_qc <- norm_beta[-unique(c(prob_out_idx, prob_noncg_idx)), -smp_out_idx]
write.table(beta_qc, paste0(dr_dir, out_file, "_beta_norm_filter.csv"), col.names = T, row.names = T)

chrXY_index <- rownames(beta_qc) %in% chrXY
pcs_qc <-  beta_qc[!chrXY_index, ]
pcs_qc <- t(na.omit(pcs_qc - rowMeans(pcs_qc)))
pc2_qc <- trlan.svd(pcs_qc, neig=10) # compute the 10 principal components
pc2_qc_u_out <- data.frame(Sample_ID = rownames(pcs_qc), 
                           PC1 = pc2_qc$u[, 1], PC2 = pc2_qc$u[, 2],
                           PC3 = pc2_qc$u[, 3], PC4 = pc2_qc$u[, 4],
                           PC5 = pc2_qc$u[, 5], PC6 = pc2_qc$u[, 6],
                           PC7 = pc2_qc$u[, 7], PC8 = pc2_qc$u[, 8],
                           PC9 = pc2_qc$u[, 9], PC10 = pc2_qc$u[, 10],
                           stringsAsFactors = F) %>%
  left_join(., pheno, by = "Sample_ID") 

pheno <- pheno %>% left_join(., pc2_qc_u_out)
write.table(pheno, paste0(dr_dir, out_file, "_Pheno.csv"), sep = ",", col.names = T, row.names = F) 


##### Export QC Plots (Raw)
grp_ctrl <- ggplot(data = ctrls_df, 
                   aes(x = 1, y = y, color = factor(CASE2), shape = COHORT)) + 
  geom_jitter(position = position_jitter(0.3), alpha = 0.7) + 
  geom_hline(aes(yintercept = threshold), color = "red", linetype = "dashed") + 
  ggtitle("BeadArray Control Reports") + 
  scale_x_continuous(limits = c(0.6, 1.4), 
                     breaks = NULL, minor_breaks = NULL) + 
  labs(x = "", y = "value", color = "Case/Control/Normal", shape = "Cohort") + 
  facet_wrap(ctrl_nm~., scales = "free_x", ncol = 3) + 
  coord_flip()


grp_sex <- ggplot(data = pred_sex, 
                  aes(x = Pred_Sex_X, y = Pred_Sex_Y, 
                      color = factor(CASE2), shape = SEX)) +
  geom_point(alpha = 0.7) +
  ggtitle("Sex Mismatches QC") + 
  scale_x_continuous(limits = c(min(pred_sex$Pred_Sex_X)-0.1, max(pred_sex$Pred_Sex_Y)+0.1)) +
  scale_y_continuous(limits = c(min(pred_sex$Pred_Sex_Y)-0.1, max(pred_sex$Pred_Sex_Y)+0.1)) +
  labs(x = "Normalized X chromosome intensities",
       y = "Normalized Y chromosome intensities",
       color = "Case/Control/Normal", shape = "Sex")


grp_geno_out <- ggplot(data = snp_outlier, 
                       aes(x = 1, y = SNP_out, 
                           color = factor(CASE2), shape = SEX)) + 
  geom_jitter(position = position_jitter(0.3), alpha = 0.7) + 
  geom_hline(yintercept = -4, color = "red", linetype = "dashed") + 
  ggtitle("Genotype calling and outliers") + 
  scale_x_continuous(limits = c(0.6, 1.4), 
                     breaks = NULL, minor_breaks = NULL) + 
  labs(y = "log odds of being outlier", x = "value", 
       color = "Case/Control/Normal", shape = "Sex") +
  coord_flip()

grp_pca <- ggplot(data = pc2_u, 
                  aes(x = PC1, y = PC2, 
                      color = factor(CASE2), shape = SEX)) +
  geom_point(alpha = 0.7) +
  ggtitle("PCA") + 
  labs(x = "PC1",
       y = "PC2",
       color = "Case/Control/Normal", shape = "Sex")

grp_lcfrac <- ggplot(data = lc_fract, 
                     aes(x = 1, y = fraction, 
                         color = factor(CASE2), shape = SEX)) + 
  geom_jitter(position = position_jitter(0.3), alpha = 0.7) + 
  ggtitle("Estimate Cell Composition") + 
  scale_x_continuous(limits = c(0.6, 1.4), 
                     breaks = NULL, minor_breaks = NULL) + 
  labs(x = "", y = "Fraction (%)", color = "Case/Control/Normal", shape = "Sex") + 
  facet_wrap(Cell~., scales = "free_x", ncol = 3) + 
  coord_flip()

pdf(file = paste0(dr_dir, out_file, "_plot_raw.pdf")) # The height of the plot in inches
print(list(grp_ctrl,
           grp_sex,
           grp_geno_out,
           grp_pca,
           grp_lcfrac))
dev.off()


##### Export QCed plots
grp_ctrl <- ggplot(data = ctrls_df %>% 
                     filter(!(Sample_ID %in% fail_smp)), 
                   aes(x = 1, y = y, color = factor(CASE2), shape = COHORT)) + 
  geom_jitter(position = position_jitter(0.3), alpha = 0.7) + 
  geom_hline(aes(yintercept = threshold), color = "red", linetype = "dashed") + 
  ggtitle("BeadArray Control Reports") + 
  scale_x_continuous(limits = c(0.6, 1.4), 
                     breaks = NULL, minor_breaks = NULL) + 
  labs(x = "", y = "value", color = "Case/Control/Normal", shape = "Cohort") + 
  facet_wrap(ctrl_nm~., scales = "free_x", ncol = 3) + 
  coord_flip()


grp_sex <- ggplot(data = pred_sex %>% 
                    filter(!(Sample_ID %in% fail_smp)), 
                  aes(x = Pred_Sex_X, y = Pred_Sex_Y, 
                      color = factor(CASE2), shape = SEX)) +
  geom_point(alpha = 0.7) + 
  scale_x_continuous(limits = c(min(pred_sex$Pred_Sex_X)-0.1, max(pred_sex$Pred_Sex_Y)+0.1)) +
  scale_y_continuous(limits = c(min(pred_sex$Pred_Sex_Y)-0.1, max(pred_sex$Pred_Sex_Y)+0.1)) +
  ggtitle("Sex Mismatches QC") + 
  labs(x = "Normalized X chromosome intensities",
       y = "Normalized Y chromosome intensities",
       color = "Case/Control/Normal", shape = "Sex")


grp_geno_out <- ggplot(data = snp_outlier %>% 
                         filter(!(Sample_ID %in% fail_smp)), 
                       aes(x = 1, y = SNP_out, 
                           color = factor(CASE2), shape = SEX)) + 
  geom_jitter(position = position_jitter(0.3), alpha = 0.7) + 
  geom_hline(yintercept = -4, color = "red", linetype = "dashed") + 
  ggtitle("Genotype calling and outliers") + 
  scale_x_continuous(limits = c(0.6, 1.4), 
                     breaks = NULL, minor_breaks = NULL) + 
  labs(y = "log odds of being outlier", x = "value", 
       color = "Case/Control/Normal", shape = "Sex") +
  coord_flip()

grp_pca <- ggplot(data = pc2_qc_u_out %>% 
                    filter(!(Sample_ID %in% fail_smp)), 
                  aes(x = PC1, y = PC2, 
                      color = factor(CASE2), shape = SEX)) +
  geom_point(alpha = 0.7) +
  ggtitle("PCA") + 
  labs(x = "PC1",
       y = "PC2",
       color = "Case/Control/Normal", shape = "Sex")

grp_lcfrac <- ggplot(data = lc_fract %>% 
                       filter(!(Sample_ID %in% fail_smp)), 
                     aes(x = 1, y = fraction, 
                         color = factor(CASE2), shape = SEX)) + 
  geom_jitter(position = position_jitter(0.3), alpha = 0.7) + 
  ggtitle("Estimate Cell Composition") + 
  scale_x_continuous(limits = c(0.6, 1.4), 
                     breaks = NULL, minor_breaks = NULL) + 
  labs(x = "", y = "Fraction (%)", color = "Case/Control/Normal", shape = "Sex") + 
  facet_wrap(Cell~., scales = "free_x", ncol = 3) + 
  coord_flip()


pdf(file = paste0(dr_dir, out_file, "_plot_out.pdf")) # The height of the plot in inches
print(list(grp_ctrl,
           grp_sex,
           grp_geno_out,
           grp_pca,
           grp_lcfrac))
dev.off()



