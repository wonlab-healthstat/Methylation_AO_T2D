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
source("/home2/sjseo/Tools/R/library/qqunif.plot.R")


#####################################
##### Set Parameters
#####################################

dr_dir <- "/home2/sjseo/Bohun/mQTL/derived_data/"
raw_snp_dir <- "/home2/sjseo/Bohun/SNParray/derived_data/"
raw_meth_dir <- "/home2/sjseo/Bohun/Methylation/derived_data/"


ncores <- detectCores()#10

#####################################
##### Summarize mQTL results
#####################################

##### Number of tests
ntest <- read.table(paste0(dr_dir, "01_01_nTest.txt"), sep = "\t", header = T, stringsAsFactors = F)
length(unique(ntest$Chr)) ; sum(ntest$nSNP) ; sum(ntest$nCpG)
sum(ntest$nTests_cis) ; sum(ntest$nTests_trans) # cis : 2,569 / trans : 407,042
with(ntest, sum(nTests_cis/nCpG))

##### SNP/CpG Info
manifest <- fread(paste0(raw_meth_dir, "manifest.csv"))
snppose <- read.table(paste0(raw_snp_dir, "IDmatched_Prune.bim"), sep = "\t", header = F, stringsAsFactors = F)
#genepos_chr <- arrange_expr %>% tibble::rownames_to_column("geneid") %>% select(geneid) %>%
#  left_join(., manifest, by = c("geneid" = "probe_id")) %>% 
#  select(geneid, chr, left = mapinfo, right = mapinfo)

##### Combine all result data
all_files <- list.files(dr_dir)

res_dat <- lapply(c("cis", "trans"), function(type) {
  files <- all_files[grep(paste0("^*(mQTL_chr)(.*)(", type, ".txt$)"), all_files)]
  dat <- lapply(files, function(f) {read.table(paste0(dr_dir, f), sep = "\t", header = T, stringsAsFactors = F)})
  dat <- do.call(rbind, dat)
  
  sig_thr <- 0.05/sum(ntest[, paste0('nTests_', type)])
  dat_dist <- dat %>% left_join(., snppose %>% select(Chr = V1, SNP = V2, pos = V4), by = "SNP") %>%
    left_join(., manifest %>% select(probe_id, mapinfo), by = c("gene" = "probe_id")) %>%
    mutate(dist = mapinfo-pos, 
           padj = p.adjust(p.value, method = "BH"),
           type = type)
})
res_dat <- do.call("rbind", res_dat)
head(res_dat)

##### Significant results of cis
res_cis <-  res_dat %>% filter(type == "cis") %>% arrange(padj) 

res_cis_rs <- lapply(unique(res_cis$Chr), function(chr) {
  tmp_dat <- res_cis %>% filter(Chr == chr)
  ref_file <-  paste0('/data/annotation/rsID_hg19_snp151/rsID_hg19_snp151_chr', chr, '.txt')
  
  print(paste("Chr", chr, "start!"))
  chr_rs <- mclapply(1:nrow(tmp_dat), function(i) {
    x <- tmp_dat[i,]
    cmd <- paste("grep -w -m 1", paste0("'", x[7], ":", x[8], "'"), ref_file, "|", "awk '{print $3}'")
    rs <- system(cmd, intern = TRUE)
    if(length(rs) == 0) rs <- NA
    return(rs)
  }, mc.cores = ncores)
  tmp_dat$rsID <- unlist(chr_rs)
  print(paste("Chr", chr, "end!"))
  
  return(tmp_dat)
})
res_cis_rs <- do.call(rbind, res_cis_rs)
write.table(res_cis_rs, paste0(dr_dir, "01_02_mQTL_cis_rsID.txt"), sep = "\t", col.names = T, row.names = F, quote = F)
#res_cis_rs <- read.table(paste0(dr_dir, "01_02_mQTL_cis_rsID.txt"), header = T,stringsAsFactors = F)
sig_res_cis <- res_cis_rs %>% filter(padj < 0.05) %>% arrange(padj)
sig_res_cis


length(unique(sig_res_cis$SNP))  #134
length(unique(sig_res_cis$gene)) #32

##### Distance vs Significance
ggplot(data = res_cis_rs, 
       aes(x = dist/1000, y = -log10(p.value), 
           color = factor(padj < 0.05, labels = c("Not sigficant", "Significant")))) +
  geom_point() + 
  geom_vline(xintercept = 0, linetype = "dashed") + 
  theme_bw() + 
  labs(x = "Distance between SNP and CpG site (kb)", 
       y = expression(-log[10](p-value)),
       color = "Significance") + 
  scale_color_manual(values = c("black", "blue")) +
  theme(legend.text = element_text(face = "bold", color = "black", size = 10))


##### Bar plot of number of SNPs
nsig_cis <- res_cis %>% group_by(gene) %>% summarize(ntol = length(SNP), nsig = sum(padj < 0.05)) %>% 
  as.data.frame() %>% mutate(notsig = ntol-nsig) 
cis_order <- nsig_cis %>% arrange(desc(ntol)) %>% mutate(order = 1:length(gene))

nsig_cis_exp <- nsig_cis %>% select(-ntol) %>%
  reshape2::melt(., id = c("gene")) %>% 
  tidyr::uncount(., weights = value) %>%
  left_join(., cis_order %>% select(gene, order), by = "gene")
nsig_cis_exp$variable <- factor(nsig_cis_exp$variable, labels = c("Significant", "Not significant"))


### Target CpG : Significant in both DM and defoliant in the methylation analysis
sig_cpg <- read.table(paste0(raw_meth_dir, "04_01_Limma_Subgroup_Table_Nested.csv"), sep = ",", header = T, stringsAsFactors = F)  
a <- sig_cpg %>% filter(Type == "defoliant_DM", significance != "Not Significant") %>%
  left_join(., manifest %>% select(probe_id, chr, mapinfo), by = c("CpG" = "probe_id"))
b <- sig_cpg %>% filter(Type == "DM", significance != "Not Significant") %>%
  left_join(., manifest %>% select(probe_id, chr, mapinfo), by = c("CpG" = "probe_id"))
target <- intersect(a$CpG, b$CpG)
cpg_target <- ifelse(unique(nsig_cis_exp$gene) %in% target, "red", "black")


ggplot(data = nsig_cis_exp, aes(x = order, fill = variable)) +
  geom_bar() + 
  annotate("text", x = cis_order$order, 
           y = cis_order$ntol+2, label =  ifelse(cis_order$nsig == 0, "", cis_order$nsig), size = 3) + 
  scale_x_continuous(breaks = unique(nsig_cis_exp$order), labels = unique(nsig_cis_exp$gene)) + 
  scale_fill_manual(values = c("blue", "grey")) +
  theme_bw() +
  theme(axis.text.x = element_text(face = "bold", color = cpg_target, 
                                   angle = 90, vjust = 0.5, hjust = 1, size =10),
        axis.text.y = element_text(face = "bold", color = "black", size = 10),
        legend.text = element_text(face = "bold", color = "black", size = 10),
        axis.line = element_line(color = "black"))+
  labs(x = "CpG", y = "Number of SNPs", fill = "")

