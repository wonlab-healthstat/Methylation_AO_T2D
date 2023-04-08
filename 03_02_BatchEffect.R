#####################################
### Study : Bohun-Methylation
### Purpose of prgram : Batch Effect Control
### Pre-required program : NA
### Output files : Start with "QCBatch"
### Programmed by : Sujin Seo
### Draft Date : 2021-10-07
### Revision : NA
### Software (version) : R (v3.6.3)
#####################################

#####################################
##### Load Required Libraries
#####################################
library(data.table)
library(dplyr)
library(lumi)
library(factoextra)
library(sva)
library(svd)
library(ggplot2)

#####################################
##### Set Parameters
#####################################
##### Derived directories
dr_dir <- "/home2/sjseo/Bohun/Methylation/derived_data/"

##### Output file
out_file <- "03_02_QCBatch"

#####################################
##### Read Data
#####################################
manifest <- fread(paste0(dr_dir, "manifest.csv"))
pheno <-  read.table(paste0(dr_dir, "03_02_QC_Pheno.csv"), sep = ",",
                     header=TRUE, stringsAsFactors=FALSE) %>%
  filter(QC_final == "PASS") %>% 
  mutate(COHORT_SEX = ifelse(COHORT == "BOHUN", 1, ifelse(SEX == "M", 2, 3)))
dat_all	<- fread(paste0(dr_dir, "03_01_QC_beta_norm_filter.csv"))
dat_all <- tibble::column_to_rownames(dat_all, "V1")
m <- beta2m(dat_all)
meth <- m[, pheno$Sample_ID]
rm(list = c("dat_all", "m"))


#####################################
##### Batch Effect Control
#####################################
##### Combat
batch <- pheno$COHORT_SEX
# modcov <- model.matrix(~CASE2, data = pheno) => cannot run
batch_cntrl1 <- ComBat(meth, batch, prior.plots = F)
save(batch_cntrl1, file = paste0(dr_dir, out_file, "_M_combat.RData"))
#load(paste0(dr_dir, out_file, "_M_combat.RData"))

##### Calculate PC
chrXY <- manifest[manifest$chr %in% c("X","Y") & 
                    manifest$probe_type!="rs", "probe_id"]
chrXY_index <- colnames(batch_cntrl1) %in% chrXY
meth_cntrl <-  t(batch_cntrl1[, !chrXY_index])
meth_cntrl_cnt <- t(na.omit(t(meth_cntrl) - rowMeans(t(meth_cntrl))))

pc2 <- trlan.svd(meth_cntrl_cnt, neig=10) # compute the first two principal components
pc2_u <- data.frame(Sample_ID = rownames(meth_cntrl), 
                    PC1 = pc2$u[, 1], PC2 = pc2$u[, 2],
                    PC3 = pc2$u[, 3], PC4 = pc2$u[, 4],
                    PC5 = pc2$u[, 5], PC6 = pc2$u[, 6],
                    PC7 = pc2$u[, 7], PC8 = pc2$u[, 8],
                    PC9 = pc2$u[, 9], PC10 = pc2$u[, 10],
                    stringsAsFactors = F)

##### Exclude Outliers
pheno_pc <- pheno %>% select(c(1:12)) %>% left_join(., pc2_u, by = "Sample_ID")  %>%
  mutate(PC_out_combat = ifelse(PC1 < 0.02 & PC2 < -0.05, "Exclude", "Include"))
write.table(pheno_pc, paste0(dr_dir, "03_02_QC_Pheno.csv"), sep = ",", col.names = T, row.names = F)

pdf(paste0(dr_dir, out_file, out_file, "_Plot_PCA.pdf"))
ggplot(data = pheno_pc %>% filter(PC_out_combat == "Include"), 
       aes(x = PC1, y = PC2, color = factor(CASE2), shape = SEX)) +
  geom_point(alpha = 0.7) +
  ggtitle("PCA") + 
  labs(x = "PC1", y = "PC2",
       color = "Case/Control/Normal", shape = "Sex")
dev.off()

res.pca <- prcomp(meth_cntrl_cnt, scale = TRUE)
pdf(paste0(dr_dir, out_file, "_Plot_PCvariables.pdf"))
fviz_pca_var(res.pca, select.var = list(contrib = 10), 
             col.var = "contrib", # Color by contributions to the PC
             #gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
)
dev.off()

