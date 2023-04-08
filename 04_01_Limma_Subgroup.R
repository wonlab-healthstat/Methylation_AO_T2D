#####################################
### Study : Bohun-Methylation
### Purpose of prgram : 
### Pre-required program : NA
### Output files : Start with "Limma_Subgroup"
### Programmed by : Sujin Seo
### Draft Date : 2021-01-20
### Revision : NA
### Software (version) : R (v4.1.1)
#####################################

#####################################
##### Load Required Libraries
#####################################
library(parallel) 
library(data.table)
library(dplyr)
library(limma)
library(ggplot2)
library(DMRcate)
library(RColorBrewer)
library(missMethyl)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(extrafont) ; loadfonts(device = "all")
library(GWASTools)
library(bacon)
#source("/home2/sjseo/Tools/R/library/qqunif.plot.R")


#####################################
##### Set Parameters
#####################################
##### Directories
dr_dir <- "/home2/sjseo/Bohun/Methylation/derived_data/"
ncores <- 30

##### Output file
out_file <- "04_01_Limma_Subgroup"


#####################################
##### Read Data
#####################################
pheno <- read.table(paste0(dr_dir, "03_02_QC_Pheno.csv"), sep = ",",
                    header=TRUE, stringsAsFactors=FALSE) %>%
  filter(QC_final == "PASS", PC_out_combat == "Include")

pheno_out <- pheno %>% select(COHORT, CASE2 = CASE, Sample_ID, SEX, DM, AGE, PC1, PC2, PC3, PC4, PC5, PC6, PC7, PC8, PC9, PC10)
dim(pheno_out)
write.table(pheno_out, paste0(dr_dir, "04_01_Phenotypes_1249.csv"), sep = ",", col.names = T, row.names = F, quote = F)
pheno <- read.table(paste0(dr_dir, "04_01_Phenotypes_1249.csv"), sep = ",", header = T)
smp <- pheno$Sample_ID
cellcomp	<- fread(paste0(dr_dir, "03_01_QC_LC.csv")) 
load(paste0(dr_dir, "03_02_QCBatch_M_combat.RData"))

manifest <- fread(paste0(dr_dir, "manifest.csv"))

## Cross-active probes
#xprobe <- read.table("~/Bohun/Methylation/data/EPIC_crossreactive.txt", header = F)
# batch_cntrl1_1 <- batch_cntrl1[!(rownames(batch_cntrl1) %in% xprobe$V1), ]
# batch_cntrl1 <- batch_cntrl1_1
# 821,509 -> 780,308

#####################################
##### Data Preparation
#####################################

##### Select only filtered subjects
dat_all_mat <- batch_cntrl1[, pheno$Sample_ID]
dat_dm_mat <- batch_cntrl1[, pheno[pheno$DM == 1, "Sample_ID"]]

##### Combine pheno and cell counts
pheno_use	<- pheno %>% 
  left_join(., cellcomp, by = c("Sample_ID")) %>%
  mutate(defoliant = ifelse(CASE2 == "Case", 1, 0))
pheno_use$CASE2 <- factor(pheno_use$CASE2, levels = c("Normal", "Control", "Case"))
pheno_use$CASE <- factor(pheno_use$CASE, levels = c("Control", "Case"))

pheno_dm <- pheno_use %>% filter(DM == 1)


#####################################
##### Limma
##### Limma Model : https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7873980/
##### Limma User Guide : http://web.mit.edu/~r/current/arch/i386_linux26/lib/R/library/limma/doc/usersguide.pdf
#####################################

##### Limma
### Defoliant among DM (Subgroup analysis)
design_mat_full <- model.matrix(~ defoliant + AGE + PC1 + PC2 + PC3 + PC4+ PC5 + PC6 + PC7 + PC8 + PC9+ PC10, 
                                data = pheno_dm) # SEX is excluded because all DM subjects are male
rownames(design_mat_full) <- pheno_dm$Sample_ID
fit_full <-  lmFit(object = dat_dm_mat, design = design_mat_full)
fit_full_eb <- eBayes(fit_full)

### Reduced (DM)
design_mat_red <- model.matrix(~ DM + AGE + SEX + PC1 + PC2 + PC3 + PC4+ PC5 + PC6 + PC7 + PC8 + PC9+ PC10, 
                               data = pheno_use)
rownames(design_mat_red) <- pheno_use$Sample_ID
fit_red <-  lmFit(object = dat_all_mat, design = design_mat_red)
fit_red_eb <- eBayes(fit_red)

### Normal vs DM w/o AO vs DM w AO
design_mat_full <- model.matrix(~ 0 + CASE2 + AGE + SEX + PC1 + PC2 + PC3 + PC4+ PC5 + PC6 + PC7 + PC8 + PC9+ PC10, 
                                data = pheno_use) 
rownames(design_mat_full) <- pheno_use$Sample_ID
colnames(design_mat_full)[1:3] <- levels(factor(pheno_use$CASE2))

fit_full <-  lmFit(object = dat_all_mat, design = design_mat_full)

contrast_mat <- makeContrasts(Case-Normal, Control-Normal, Case-Control,
                              levels = design_mat_full)
fit_full_ctrst <- contrasts.fit(fit_full, contrasts = contrast_mat)
fit_full_ctrst_eb <- eBayes(fit_full_ctrst)

save(list= c("fit_full_eb", "fit_red_eb", "fit_full_ctrst_eb"), 
     file = paste0(dr_dir, out_file, "_Table_Nested.RData"))
#load(file = paste0(dr_dir, out_file, "_Table_Nested.RData"))


##### Summarize Results

### Different in any groups
head(topTableF(fit_full_ctrst_eb))
fit_any_res <- topTableF(fit_full_ctrst_eb, number=Inf, p.value=1, adjust.method = "BH") %>%
  tibble::rownames_to_column("CpG")

### Different in Case-Normal
fit_dm_res <- topTable(fit_full_ctrst_eb, coef = 1, number=Inf, p.value=1, adjust.method = "BH") %>%
  tibble::rownames_to_column("CpG") 

### Different in Control-Normal
fit_dmwo_res <- topTable(fit_full_ctrst_eb, coef = 2, number=Inf, p.value=1, adjust.method = "BH") %>%
  tibble::rownames_to_column("CpG") 

### Different in Case-Control
fit_def_res  <- topTable(fit_full_ctrst_eb, coef = 3, number=Inf, p.value=1, adjust.method = "BH") %>%
  tibble::rownames_to_column("CpG") 


### Different in DM - Normal
result_dm <- limma::topTable(fit_red_eb, coef = "DM", number=Inf, p.value=1) %>%
  tibble::rownames_to_column(., var = "CpG") %>%
  mutate(Type = "DM")


### Different in Case-Control among DM
result_dm_def <- limma::topTable(fit_full_eb, coef = "defoliant", number=Inf, p.value=1) %>%
  tibble::rownames_to_column(., var = "CpG") %>%
  mutate(Type = "defoliant_DM")
         
png(paste0(dr_dir, out_file, "_Plot_QQplot.png"))
par(mfrow = c(1,2))
qqPlot(result_dm$P.Value)
qqPlot(result_dm_def$P.Value)
par(mfrow = c(1,1))
dev.off()

# ########## For genomic control
# ### Different in DM - Normal
# result_dm <- limma::topTable(fit_red_eb, coef = "DM", number=Inf, p.value=1) %>%
#   tibble::rownames_to_column(., var = "CpG") %>%
#   mutate(Type = "DM", 
#          z = qnorm(P.Value/2, lower.tail = F)*sign(t),
#          chistat = z^2,
#          lambda = median(chistat)/0.455,
#          gc_stat = chistat/lambda,
#          gc_pval = 1-pchisq(gc_stat, 1))
# result_dm_gc <- bacon(result_dm$z)
# result_dm <- result_dm %>% 
#   mutate(.,
#          bc_z = tstat(result_dm_gc), 
#          bc_P.Val =  pval(result_dm_gc),
#          bc_val_adj.P.Val = p.adjust(bc_P.Val, method = "BH"))
# 
# 
# ### Different in Case-Control among DM
# result_dm_def <- limma::topTable(fit_full_eb, coef = "defoliant", number=Inf, p.value=1) %>%
#   tibble::rownames_to_column(., var = "CpG") %>%
#   mutate(Type = "defoliant_DM",
#          z = qnorm(P.Value/2, lower.tail = F)*sign(t),
#          chistat = z^2,
#          lambda = median(chistat)/0.455,
#          gc_stat = chistat/lambda,
#          gc_pval = 1-pchisq(gc_stat, 1))  
# result_dm_def_gc <- bacon(result_dm_def$z)
# result_dm_def <- result_dm_def %>% 
#   mutate(.,
#          bc_z = tstat(result_dm_def_gc),
#          bc_P.Val =  pval(result_dm_def_gc),
#          bc_val_adj.P.Val = p.adjust(bc_P.Val, method = "BH"))
# 
# par(mfrow = c(1,2))
# qqPlot(result_dm$bc_P.Val)
# qqPlot(result_dm_def$bc_P.Val)
# par(mfrow = c(1,1))
# 
# par(mfrow = c(1,2))
# qqPlot(result_dm$gc_pval)
# qqPlot(result_dm_def$gc_pval)
# par(mfrow = c(1,1))


# x2 <- rchisq(100,1,.1)
# p <- pchisq(x2,1,lower.tail=FALSE)
# r <- gcontrol2(p)
# print(r$lambda)



##### Volcano Plot
fc_thr <- 0.1

result <- rbind(result_dm, result_dm_def) %>%
  mutate(significance = ifelse(logFC > fc_thr & adj.P.Val <= 0.05, "Hyper", 
                               ifelse(logFC < -fc_thr & adj.P.Val <= 0.05 , "Hypo", "Not Significant"))) 
result$significance <- factor(result$significance, levels = c("Not Significant", "Hypo", "Hyper"))
write.table(result, paste0(dr_dir, out_file, "_Table_Nested.csv"), sep = ",", row.names = F, col.names = T)  
result <- read.csv(paste0(dr_dir, out_file, "_Table_Nested.csv"), header= T)

#result <- read.table(paste0(dr_dir, out_file, "_Table_Nested.csv"), sep = ",", header = T, stringsAsFactors = F)
result%>% mutate(sig_pval = ifelse(adj.P.Val <= 0.05, "sig", "not-sig"),
                 sig_fc = ifelse(abs(logFC) > fc_thr, "sig", "not_sig")) %>%
  dplyr::group_by(Type, sig_pval, sig_fc) %>% dplyr::summarize(N = length(CpG))


type <- c("DM", "defoliant_DM")
max(abs(result$logFC))            # To scale x-axis
max(-log10(result$adj.P.Val))     # To scale y-axis
for (x in type) {
  png(paste0(dr_dir, out_file, "_Plot_Volcano_", x, ".png"), width = 300, height = 250)
  plt <- ggplot(result %>% filter(Type == x), aes(logFC, -log10(adj.P.Val), col = significance)) +
    geom_point() +
    geom_vline(xintercept = c(-fc_thr, fc_thr), linetype =2) +
    scale_x_continuous(limits = c(-max(abs(result$logFC))-0.5, max(abs(result$logFC))+0.5)) +
    scale_y_continuous(limits = c(0, max(-log10(result$adj.P.Val)) +1)) +
    scale_color_manual(values = c("black", "blue", "red")) +
    labs(x = expression(log[2](FC)), y = expression(-log[10](p-value))) +
    theme_bw() +
    theme(legend.position = "none", text=element_text(size=13,  family="Arial"),
          axis.line.y.left = element_line(color = "black"),
          axis.line.y.right = element_line(color = "black"),
          axis.line.x.top = element_line(color = "black"),
          axis.line.x.bottom = element_line(color = "black"))
  print(plt)
  dev.off()
}


### Target CpG
a <- result %>% filter(Type == "DM", abs(logFC) > fc_thr, adj.P.Val < 0.05)
b <- result %>% filter(Type == "defoliant_DM", abs(logFC) > fc_thr, adj.P.Val < 0.05)
target <- intersect(a$CpG, b$CpG)

### Plot of Changes of significant CpGs
dm <- fit_dmwo_res[fit_dmwo_res$CpG %in% target, ] %>% #tibble::rownames_to_column("CpG") %>% 
  mutate(Type = "DMwoAO")
def_dm2 <- fit_def_res[fit_def_res$CpG %in% target, ] %>% #tibble::rownames_to_column("CpG") %>% 
  mutate(Type = "AO")
def_dm <- fit_dm_res[fit_dm_res$CpG %in% target, ] %>% #tibble::rownames_to_column("CpG") %>% 
  mutate(Type = "DMwAO")
both <- do.call(rbind, list(dm, def_dm2, def_dm))
both$Type <- factor(both$Type, levels = c("DMwoAO", "AO", "DMwAO"),  
                    labels = c("AO-unexposed T2D\nvs Normal", 
                               "AO-exposed T2D\nvs AO-unexposed T2D",
                               "AO-exposed T2D\nvs Normal"))
arrow <- dm %>% select(CpG, FromX = logFC, FromY = P.Value) %>%
  left_join(., def_dm2 %>% select(CpG, ToX = logFC, ToY = P.Value), by = "CpG")
arrow2 <- def_dm2 %>% select(CpG, FromX = logFC, FromY = P.Value) %>%
  left_join(., def_dm %>% select(CpG, ToX = logFC, ToY = P.Value), by = "CpG")

text_axis_x <- def_dm$logFC+0.05
text_axis_y <- -log10(def_dm$P.Value)+2
  
png(paste0(dr_dir, out_file, "_Plot_SigChanges.png"), width = 800, height = 400)
ggplot(data = both, 
       aes(x = logFC, y = -log10(P.Value), color = Type, shape = Type)) +
  geom_point(size = 2.5) +
  geom_segment(data = arrow, aes(x = FromX, y = -log10(FromY), xend = ToX, yend = -log10(ToY)),
               #arrow = arrow(length = unit(0.2, "cm")),
               color = "darkgrey",
               inherit.aes = F) +
  geom_segment(data = arrow2, aes(x = FromX, y = -log10(FromY), xend = ToX, yend = -log10(ToY)),
               arrow = arrow(length = unit(0.2, "cm")),
               color = "darkgrey",
               inherit.aes = F) +
  ggrepel::geom_text_repel(data = both %>% filter(Type == "AO-exposed T2D\nvs Normal"),
                           aes(x = logFC, y = -log10(P.Value), label = CpG), 
                           nudge_x = 0, nudge_y = 1.5, segment.alpha = 0, 
                           size = 3.5, hjust = 0, family="Arial", inherit.aes = F) + 
  theme_bw() +
  scale_shape_manual("", values = c(1, 16, 16)) +
  scale_color_manual("", values = c("black", grey.colors(2, start = 0.4, end = 0))) +
  scale_x_continuous(limits = c(-0.5, 0.55)) + 
  labs(x = expression(log[2](FC)), y = expression(-log[10](p-value)), color = "") +
  theme(text=element_text(size=13,  family="Arial", color = "black"),
        legend.key.height = unit(1.5, "cm"),
        axis.line = element_line(colour = 'black', size = 0.2))

dev.off()




### DM vs Normal / DM wo AO vs DM w AO
target <- setdiff(a$CpG, b$CpG)

### Plot of Changes of significant CpGs
dm <- result_dm[result_dm$CpG %in% target, ] %>% #tibble::rownames_to_column("CpG") %>% 
  mutate(Type = "DM")
def_dm <- fit_def_res[fit_def_res$CpG %in% target, ] %>% #tibble::rownames_to_column("CpG") %>% 
  mutate(Type = "DMwAO")
both <- rbind(dm, def_dm)
both$Type <- factor(both$Type, levels = c("DM", "DMwAO"),  
                    labels = c("DM vs Normal", 
                               "AO-exposed T2D\nvs AO-unexposed T2D"))
arrow <- dm %>% select(CpG, FromX = logFC, FromY = P.Value) %>%
  left_join(., def_dm %>% select(CpG, ToX = logFC, ToY = P.Value), by = "CpG")

png(paste0(dr_dir, out_file, "_Plot_NotSigChanges.png"), width = 700, height = 400)
ggplot(data = both, 
       aes(x = logFC, y = -log10(P.Value), color = Type)) +
  geom_vline(xintercept = 0, linetype = 'dashed') + 
  geom_segment(data = arrow, aes(x = FromX, y = -log10(FromY), xend = ToX, yend = -log10(ToY)),
               arrow = arrow(length = unit(0.2, "cm")),
               color = "darkgrey",
               inherit.aes = F) +
  scale_color_manual(values = c("darkgrey", "black")) +
  geom_point(size = 2.5) +
  theme_bw() +
  labs(x = expression(log[2](FC)), y = expression(-log[10](p-value)), color = "")+
  theme(text=element_text(size=13,  family="Arial"))
dev.off()


### Number of CpG site which intersects each other
result %>% group_by(Type, significance) %>% summarize(N = length(CpG))

sig_res1 <- result %>% filter(Type == "DM") %>% 
  dplyr::select(CpG, DM_fc = logFC, DM_sig = significance)
sig_res2 <- result %>% filter(Type == "defoliant_DM") %>% 
  dplyr::select(CpG, DEF_fc = logFC, DEF_sig = significance)
sig_res <- sig_res1 %>% left_join(., sig_res2, by = "CpG") %>%
  mutate(sig_change = paste(DM_sig, "->", DEF_sig)) %>%
  left_join(., manifest %>% select(probe_id, chr, addressU), by = c("CpG" = "probe_id"))

with(sig_res, table(DM_sig, DEF_sig))
sig_res %>% filter(CpG %in% target)


#####################################
##### Gene Set Analysis 
##### Gene ontology terms for biological process
#####################################
sig_cpg <- target
gst_go <- gometh(sig.cpg = sig_cpg, array.type = "EPIC",
                 all.cpg = sig_res$CpG, 
                 collection = "GO", plot.bias = FALSE)
gst_kegg <- gometh(sig.cpg = sig_cpg, array.type = "EPIC",
                   all.cpg = sig_res$CpG, 
                   collection = "KEGG")
gst_sig <- list(gst_go %>% filter(FDR < 0.05), 
                gst_kegg %>% filter(FDR < 0.05))
lapply(gst_sig, dim) ## None are found



#####################################
##### DMR
#####################################

##### Remove X and Y chromosome 
noSNPs <- rmSNPandCH(as.matrix(dat_all_mat), 
                     dist=2, mafcut=0.05, rmXY = TRUE)

### Defoliant|DM
annot_def_dm <- cpg.annotate(datatype = "array", 
                             object = as.matrix(dat_dm_mat), what="M", 
                             arraytype="EPIC",
                             design = design_mat_full, 
                             coef = "defoliant1")

dmr_def_dm <- dmrcate(annot_def_dm)

dmr_def_dm_ranges <- extractRanges(dmr_def_dm, genome="hg19")
dmr_def_dm_df <- data.frame(dmr_def_dm_ranges, stringsAsFactors=FALSE) %>% 
  mutate(Type = "def_DM")

### DM
annot_dm <- cpg.annotate(datatype = "array", 
                         object = as.matrix(dat_all_mat), what="M", 
                         arraytype="EPIC",
                         design = design_mat_red , 
                         coef = "DM1")

dmr_dm <- dmrcate(annot_dm)
dmr_dm_ranges <- extractRanges(dmr_dm, genome="hg19")
dmr_dm_df <- data.frame(dmr_dm_ranges, stringsAsFactors=FALSE) %>%
  mutate(Type = "DM")

#png(paste0(dr_dir, out_file, "_DMR.png"), width=2000, height=2000)
#DMR.plot(ranges = res_ranges, 
#         dmr = which.min(res2$min_smoothed_fdr), 
#         CpGs = as.matrix(dat_all_mat), what="M",
#         arraytype="EPIC", genome="hg19", samps=1:6, 
#         phen.col = rep("forestgreen", nrow(dat_all_mat)))
#dev.off() 
