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
library(ivtools)
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

#ncores <- detectCores()
ncores <- 20


#########################################################
##### Mendelian Randomization
#########################################################
find_instruments <- function(formula1, formula2) {
  f1_covariates <- all.vars(formula1)[-1]
  instrumented <- f1_covariates[f1_covariates %in% all.vars(formula2)]
  instruments <- all.vars(formula2)[-1]
  out <- list(instruments = instruments,
              instrumented = instrumented)
  return(out)
}

iv.glm <- function(model_formula,
                   instrument_formula = NULL, data = NULL,
                   family = binomial, link = 'logit', ...) {
  if(is.null(instrument_formula)) {
    out <- glm(model_formula, data = data, family = family(link = link))
  } else {
    if(class(model_formula) != "formula" | class(instrument_formula) != "formula") {
      stop("both model_formula and instrument_formula must be of class formula")
    }
    # Instrument info
    instrument_data <- find_instruments(model_formula, instrument_formula)
    instruments <- instrument_data$instruments
    instrumented <- instrument_data$instrumented
    n_instruments <- length(instruments)
    if(length(instrumented) > 1) {
      stop("You may only instrument one variable at a time.")
    }
    
    # Estimate model
    stage_one <- lm(instrument_formula, data = model.frame(instrument_formula, data = data))
    model_data <- model.frame(model_formula, data = data)
    model_data[, which(names(model_data) == instrumented)] <- stage_one$fitted.values
    stage_two <- glm(formula = model_formula, data = model_data,
                     family = family(link = link))
    
    # Diagnostics
    cor_w_var <- vector('numeric', length = n_instruments)
    cor_w_error <- vector('numeric', length = n_instruments)
    for(i in seq_along(1:n_instruments)) {
      cor_w_var[i] <- cor(stage_one$model[, 1], stage_one$model[, (i+1)])
      cor_w_error[i] <- cor(stage_one$model[, (i+1)], stage_two$residuals)
      cat(paste0("correlation between ", names(stage_one$model)[1], " and ", names(stage_one$model)[(i+1)], ": ",
                 round(cor_w_var[i], 3),  "\n"))
      cat(paste0("correlation between ", names(stage_one$model)[(i+1)], " and residuals: ",
                 round(cor_w_error[i], 3),  "\n"))
    }
    
    # Return results and diagnostics
    out <- list(
      exclusion_restriction = round(cor_w_error, 3),
      instrument_validity = round(cor_w_var, 3),
      instruments = instruments,
      instrumented = instrumented,
      stage_one = stage_one,
      fit = stage_two
    )
    class(out) <- "ivm"
  }
  return(out)
}



#########################################################
##### Mendelian Randomization (Single)
#########################################################

##### Match sample ID with barcode ID
meth_smp1 <- read.table(file.path(raw_meth_dir2, "Sample.Table.txt"), sep = "\t", header = T, stringsAsFactors = F)
meth_smp2 <- read.table(file.path(raw_meth_dir3, "Sample.Table 1.txt"), sep = "\t", header = T, stringsAsFactors = F)
meth_smp3 <- read.table(file.path(raw_meth_dir3, "Sample.Table 2.txt"), sep = "\t", header = T, stringsAsFactors = F)
meth_smp <- rbind(meth_smp1, meth_smp2, meth_smp3) %>%
  mutate(Barcode = paste(Sentrix.Barcode, Sample.Section, sep = "_")) %>% 
  select(Barcode, Sample.ID)

pheno <- read.table(paste0(raw_meth_dir, "04_01_Phenotypes_1249.csv"), sep = ",",
                    header=TRUE, stringsAsFactors=FALSE) %>%
  #filter(QC_final == "PASS", PC_out_combat == "Include") %>%
  left_join(., meth_smp, by = c("Sample_ID" = "Barcode")) %>%
  mutate(Sample_ID = ifelse(is.na(Sample.ID), Sample_ID, Sample.ID), 
         defoliant = ifelse(CASE2 == "Case", 1, 0), 
         CASE = ifelse(defoliant == 1, "Case", ifelse(DM == 1, "Control", "Normal"))) %>%
  select(-Sample.ID)

load(paste0(dr_dir, "mQTL_methyl_dat.RData"))
meth <- batch_cntrl1 ; rm(batch_cntrl1)

res_cis_rs <- read.table(paste0(dr_dir, "01_02_mQTL_cis_rsID.txt"), header = T,stringsAsFactors = F)
snp_list <- res_cis_rs %>% filter(padj < 0.05) %>% arrange(padj)


#final_target <- snp_list %>% filter(gene %in% target, padj < 0.05)
#sem <- lapply(1:nrow(snp_list), function(i) {
sem <- mclapply(1:nrow(snp_list), function(i) {
  print(i)
  
  target_cpg <- snp_list[i,"gene"]
  target_snp <- snp_list[i,"SNP"]
  tmp_rsID <- snp_list[i,"rsID"]
  tmp_chr <- snp_list[i, 'Chr']
  tmp_pos <- snp_list[i, 'pos']
  tmp_mapinfo <- snp_list[i, 'mapinfo']
  
  ### Causation
  tmp_snp_dat <- t(read.table(paste0(dr_dir, "01_01_SNP/SNP_chr", tmp_chr, ".txt"))[target_snp, ])
  tmp_meth_dat <- t(meth[target_cpg, rownames(tmp_snp_dat)])
  tmp_dat <- data.frame(SNP = tmp_snp_dat, 
                        METH = tmp_meth_dat,
                        stringsAsFactors = F) %>%
    tibble::rownames_to_column("Sample_ID") %>%
    left_join(., pheno, by = "Sample_ID") %>%
    mutate(CASE2_NUM = as.numeric(factor(CASE2, levels = c("Control", "Case")))-1)
  tmp_dat <- tmp_dat[complete.cases(tmp_dat),]
  colnames(tmp_dat)[2:3] <- c("SNP", "Meth")
 
  tmp_dat$CASE2 <- factor(tmp_dat$CASE2, levels = c("Control", "Case"))
  if(all(table(tmp_dat$CASE2)) > 0) {
    fitX.LZ <- glm(formula=Meth~SNP, data=tmp_dat)
    fitY.LX <- glm(formula=CASE2~Meth+AGE+SEX, 
                   family="binomial", data=tmp_dat)
    fitIV <- ivglm(estmethod="ts", fitX.LZ=fitX.LZ, fitY.LX=fitY.LX, data=tmp_dat, ctrl=TRUE) 
    
    ive_coef0 <- (coef(summary(fitX.LZ)))[2,]
    ive_coef <- (coef(summary(fitIV)) %>% data.frame())[5,]
    ive_coef_res <- data.frame(CHR = tmp_chr, SNP = tmp_rsID, SNP_lab = target_snp, POS = tmp_pos,
                               Meth = target_cpg, MAPINFO = tmp_mapinfo,
                               SNP_Meth_Est = ive_coef0[1], SNP_Meth_SE = ive_coef0[2], SNP_Meth_Stat = ive_coef0[3], SNP_Meth_Pval = ive_coef0[4],
                               Meth_Pheno_Est = ive_coef[, 1], Meth_Pheno_SE = ive_coef[, 2], Meth_Pheno_Stat = ive_coef[, 3], Meth_Pheno_Pval = ive_coef[, 4])
    
    return(ive_coef_res)
  }
  }, mc.cores = ncores)

#})

# https://cran.r-project.org/web/packages/ivreg/vignettes/ivreg.html
# https://github.com/alicerosecarter/MediationMR/blob/master/example_R.R


sem_res <- do.call(rbind, lapply(sem, function(x) x)) %>%
  mutate(SNP_Meth_Adj_Pval = p.adjust(SNP_Meth_Pval), 
         Meth_Pheno_Adj_Pval = p.adjust(Meth_Pheno_Pval)) %>% 
  arrange(Meth_Pheno_Pval)


ao_cpg <- c("cg07553761", "cg20075319", "cg21757266", "cg05203217", "cg20102280", "cg26081717", "cg21878650")
unique(sem_res %>% filter(Meth_Pheno_Adj_Pval < 0.05) %>% select(Meth))
sem_res %>% filter(Meth %in% ao_cpg)
write.table(sem_res, paste0(dr_dir, "02_01_MR_single.txt"), sep = "\t", 
            col.names = T, row.names = F, quote = F)


#########################################################
##### Boxplot (After 02_02_MR_reverse.R)
#########################################################
#box_dat <- logit_res2 %>% filter(Causal == "Reactive")

sem_res <- read.table(paste0(dr_dir, "02_01_MR_single.txt"), sep = "\t", header = T, stringsAsFactors = F)
sem_reverse <- read.table(paste0(dr_dir, "02_02_MR_reverse_single_0927.txt"), sep = "\t", header = T, stringsAsFactors = F)
sem_reverse %>% filter(Pheno_Meth_Adj_Pval < 0.05)
sem_reverse %>% filter(Meth %in% ao_cpg)

box_dat <- sem_res %>% filter(Meth_Pheno_Pval < 0.05) %>% 
  #sem_reverse %>% filter(Pheno_Meth_Pval < 0.05, SNP_lab %in% gwas_snp) %>% 
  mutate(SNP = as.character(SNP),
         SNP_lab = as.character(SNP_lab))


sig_chr <- unique(box_dat$CHR)
# sig_chr <- c(1, 9, 13)
bxplt <- lapply(sig_chr, function(chr) {
  
  target_cpg <- unique(box_dat[box_dat$CHR == chr,"Meth"])
  target_snp <- box_dat[box_dat$CHR == chr,"SNP"]
  snp_nm <- box_dat[box_dat$CHR == chr, "SNP_lab"]
  
  #tmp_chr <- box_dat[i, "CHR"]
  tmp_snp_dat <- t(read.table(paste0(dr_dir, "01_01_SNP/SNP_chr", chr, ".txt"), stringsAsFactors = F)[snp_nm, ]) %>%
    reshape2::melt(.) %>% dplyr::rename(Sample_ID = Var1, SNP_lab = Var2, allele = value) %>%
    left_join(., box_dat %>% select(SNP, SNP_lab), by = "SNP_lab")
  tmp_snp_dat$Sample_ID <- as.character(tmp_snp_dat$Sample_ID)
  tmp_snp_dat$SNP <- as.character(tmp_snp_dat$SNP)
  tmp_meth_dat <- t(meth[target_cpg, tmp_snp_dat$Sample_ID]) %>% 
    reshape2::melt(.) %>% dplyr::rename(Sample_ID = Var1, Meth = Var2, Beta = value)
  tmp_meth_dat$Beta <- 2^tmp_meth_dat$Beta / (2^tmp_meth_dat$Beta + 1)
  tmp_dat <- tmp_snp_dat %>% left_join(., tmp_meth_dat, by = "Sample_ID") %>%
    left_join(., pheno %>% select(Sample_ID, CASE), by = "Sample_ID")
  
  clr <- RColorBrewer::brewer.pal(3, 'Blues')
  plt <- ggplot(data = na.omit(tmp_dat), 
                aes(x = factor(allele), y = Beta,
                    fill = factor(CASE, levels = c("Case", "Control", "Normal"),
                                  labels = c("AO-exposed T2D", "AO-unexposed T2D", "Normal")))) +
    geom_boxplot() +
    labs(x = "Number of minor alleles", y = "Methylation (Beta)\n", fill = "Case") +
    scale_fill_manual(values = c(clr[3], clr[2], clr[1])) +
    theme_bw() +
    facet_grid(Meth~SNP)
  
  return(plt)
})

legend <-get_legend(bxplt[[1]] + theme(legend.position = "bottom"))

all_plt <- ggdraw() +
  draw_plot(bxplt[[1]] +  theme(legend.position = "none"), x = 0, y = 0.5, width = 1, height = 0.5) +
  draw_plot(bxplt[[2]] +  theme(legend.position = "none"), x = 0, y = 0, width = 0.7, height = .5) +
  draw_plot(bxplt[[3]] +  theme(legend.position = "none"), x = .7, y = 0, width = 0.3, height = .5)
plot_grid(all_plt, legend, ncol = 1, rel_heights = c(1, .1))  

