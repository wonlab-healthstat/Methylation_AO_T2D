# PT and CT
library(data.table); library(stringr); library(parallel)
# LDpred
library(data.table); library(magrittr); library(bigsnpr); library(parallel)
options(bigstatsr.check.parallel.blas=FALSE); options(default.nproc.blas=NULL)
# LASSOsum
library(data.table); library(lassosum); library(fdrtool); library(parallel)
# PRScs
library(data.table); library(parallel)

# Function definition
uniPRS.fun = function(name.idir,name.ss,name.test,name.odir,method){
  
  print(method)
  
  # Specify the test data and read summary statistics 
  test.bfile = paste0(name.idir,name.test)
  ss = fread(paste0(name.idir,name.ss,".sumstats"),data.table=F)
  fam.df = fread(paste0(test.bfile,".fam"),data.table=F)
  ncore = detectCores()
  
  if(method=="PT"|method=="CT"){
    
    cutoff = 10^(-5)
    
    if(method=="PT"){
      
      # Pruning and Extract the SNPs 
      ss.df = ss[ss$P<cutoff,]
      system(paste("plink --bfile",test.bfile,"--indep-pairwise 50 5 0.5 --out",test.bfile,"--threads",ncore))
      test.pr.bfile = paste0(test.bfile,"_pr0.5")
      system(paste("plink --bfile",test.bfile,"--extract",paste0(test.bfile,".prune.in"),"--make-bed --out",test.pr.bfile,"--threads",ncore))
      
      bim = fread(paste0(test.pr.bfile,".bim"),data.table=F)
      bim.df = bim[bim$V2%in%ss.df$SNP,]
      ss.df = ss.df[ss.df$SNP%in%bim.df$V2,]
      name.snp = paste0(test.bfile,"_pruned_SNP.txt")
      write.table(data.frame(bim.df$V2),file=name.snp,row.names=F,col.names=F,quote=F)
      
      
      test.bfile.pr = paste0(test.bfile,"_pr0.5")
      system(paste("plink --bfile",test.pr.bfile,"--extract",name.snp,"--recodeA --freq --out",test.bfile.pr,"--threads",min(80,ncore)))
      raw = fread(paste0(test.bfile.pr,".raw"),data.table=F)
      maf.df = fread(paste0(test.bfile.pr,".frq"),data.table=F)
      
    }else if(method=="CT"){
      
      # Clumping and Extract the SNPs
      system(paste("plink --bfile",test.bfile,"--clump",paste0(name.idir,name.ss,".sumstats"),"--clump-p1",cutoff,"--clump-p2",cutoff,"--clump-r2 0.5 --out",test.bfile,"--threads",ncore))
      clp.df = fread(paste0(test.bfile,".clumped"),data.table=F)
      ss.df = ss[ss$SNP%in%clp.df$SNP,]
      bim = fread(paste0(test.bfile,".bim"),data.table=F)
      ss.df = ss.df[ss.df$SNP%in%bim$V2,]
      bim.df = bim[bim$V2%in%ss.df$SNP,]
      match.ss.vec = match(bim.df$V2,ss.df$SNP)
      ss.df = ss.df[match.ss.vec,]
      name.snp = paste0(test.bfile,"_clumped_SNP.txt")
      write.table(data.frame(SNP=bim.df$V2),file=name.snp,row.names=F,col.names=F,quote=F)
      test.bfile.clp = paste0(test.bfile,"_clp0.5")
      system(paste("plink --bfile",test.bfile,"--extract",name.snp,"--recodeA --freq --out",test.bfile.clp,"--threads",ncore))
      raw = fread(paste0(test.bfile.clp,".raw"),data.table=F)
      maf.df = fread(paste0(test.bfile.clp,".frq"),data.table=F)
      
    }
    
    # Fill the NAs in .raw
    id.df = raw[,c("FID","IID")]; raw.df = raw[,-(1:6)]
    min.vec = unlist(lapply(strsplit(colnames(raw.df),split="_"),function(X){ X[2] }))
    raw.vec = unlist(lapply(strsplit(colnames(raw.df),split="_"),function(X){ X[1] })); colnames(raw.df) = raw.vec
    if(sum(is.na(raw.df))!=0){
      imp.mat = 2*t(matrix(rep(maf.df$MAF,nrow(id.df)),ncol=nrow(id.df)))
      raw.df[is.na(raw.df)] = imp.mat[is.na(raw.df)]
    }
    
    # Match alleles in .raw to sumstats : (ss$REF-bim$V6)&(ss$ALT-bim$V5) should be equal
    if(sum(ss.df$ALT!=min.vec)!=0){
      rm.vec = c()
      for(pos.iter in 1:ncol(raw.df)){
        raw.min = min.vec[pos.iter]; ss.min = ss.df[pos.iter,"ALT"]
        raw.maj = bim.df$V6[pos.iter]; ss.maj = ss.df[pos.iter,"REF"]
        
        if((raw.min==ss.min)+(raw.maj==ss.maj)==2){
          raw.df[,pos.iter] = raw.df[,pos.iter]
        }else if((raw.min==ss.min)+(raw.maj==ss.maj)==1){
          rm.vec = c(rm.vec,colnames(raw.df)[pos.iter])
        }else if((raw.min==ss.min)+(raw.maj==ss.maj)==0){
          if((raw.min==ss.maj)*(raw.maj==ss.min)==1){ raw.df[,pos.iter] = 2-raw.df[,pos.iter] }else{ rm.vec = c(rm.vec,colnames(raw.df)[pos.iter]) }
        }
        
      }
      if(length(rm.vec)!=0){ raw.df = raw.df[,!(colnames(raw.df)%in%rm.vec)]; ss.df = ss.df[!(ss.df$SNP%in%rm.vec),] }
    }
    
    # Calculate PRS 
    beta.mat = t(matrix(rep(ss.df$BETA,nrow(raw.df)),ncol=nrow(raw.df)))
    prs.vec = rowSums(as.matrix(raw.df)*beta.mat)
    prs.df = data.frame(FID=fam.df$V1,IID=fam.df$V2,PRS=prs.vec)
    write.table(prs.df,file=paste0(name.odir,name.test,".",method,".prs"),row.names=F,quote=F)
    
  }
  
  if(method=="LDpred"){
    
    if(!dir.exists(paste0(name.idir,"LDpred-tmp"))){ system(paste("mkdir",paste0(name.idir,"LDpred-tmp"))) }
    # Zip the summary statistics & Amend the column names
    if(!file.exists(paste0(name.idir,"LDpred-tmp/",name.ss,".sumstats.gz"))){
      system(paste("cp",paste0(name.idir,name.ss,".sumstats"),
                   paste0(name.idir,"LDpred-tmp/",name.ss,".sumstats")))
      system(paste("gzip",paste0(name.idir,"LDpred-tmp/",name.ss,".sumstats")))
    }
    ss = bigreadr::fread2(paste0(name.idir,"LDpred-tmp/",name.ss,".sumstats.gz"))
    ss = ss[,c("CHR","BP","SNP","ALT","REF","N","SE","P","BETA")]
    names(ss) = c("chr","pos","rsid","a1","a0","n_eff","beta_se","p","beta")
    
    #ss = ss[,c("CHR","BP","SNP","ALT","REF","N","P","BETA")]
    #names(ss) = c("chr","pos","rsid","a1","a0","n_eff","p","beta")
    
    # Obtain HapMap3 SNPs & Filter out hapmap SNPs
    info = readRDS(url("https://ndownloader.figshare.com/files/25503788"))
    ss.df = ss[ss$rsid%in%info$rsid,]
    
    # Calculate the LD matrix
    tmp = tempfile(tmpdir=paste0(name.idir,"LDpred-tmp")) # Open a temporary file
    on.exit(file.remove(paste0(tmp,".sbk")),add=TRUE)
    corr = NULL; ld = NULL # initialize variables for storing the LD score and LD matrix
    
    fam.order = NULL # order of samples in the bed file
    if(!file.exists(paste0(test.bfile,".bk"))){ snp_readBed(paste0(test.bfile,".bed")) }
    obj.bigSNP = snp_attach(paste0(test.bfile,".rds")) # attach the genotype object
    
    map.df = obj.bigSNP$map[-3] # extract the SNP information from the genotype
    names(map.df) = c("chr", "rsid", "pos", "a1", "a0")
    info_snp = snp_match(ss.df,map.df,match.min.prop=0.5) # perform SNP matching
    genotype = obj.bigSNP$genotypes
    CHR = map.df$chr; POS = map.df$pos
    
    POS2 = snp_asGeneticPos(CHR,POS,dir=paste0(name.idir,"LDpred-tmp")) # get the CM information from 1000 Genome
    for(chr in 1:22){ # calculate LD
      # Extract SNPs that are included in the chromosome
      ind.chr = which(info_snp$chr == chr)
      ind.chr2 = info_snp$"_NUM_ID_"[ind.chr]
      corr0 = snp_cor(genotype,ind.col=ind.chr2,ncores=min(80,ncore),infos.pos=POS2[ind.chr2],size=3/1000)
      if(chr==1){
        ld = Matrix::colSums(corr0^2)
        corr = as_SFBM(corr0, tmp)
      }else{
        ld = c(ld, Matrix::colSums(corr0^2))
        corr$add_columns(corr0, nrow(corr))
      }
    }
    
    fam.order = as.data.table(obj.bigSNP$fam)
    setnames(fam.order,c("family.ID","sample.ID"),c("FID","IID"))
    
    # Perform LD score regression
    df_beta = info_snp[,c("beta","beta_se","n_eff","_NUM_ID_")]
    ldsc = snp_ldsc(ld,length(ld),chi2=(df_beta$beta/df_beta$beta_se)^2,sample_size=df_beta$n_eff,blocks = NULL)
    h2_est = ldsc[["h2"]]
    
    # Obtain LDpred adjusted beta & PRS 
    if(is.null(obj.bigSNP)){ obj.bigSNP = snp_attach(paste0(test.bfile,".rds")) }
    genotype = obj.bigSNP$genotypes
    ind.test = 1:nrow(genotype)
    
    for(opt.iter in c("inf","grid","auto")){ 
      
      print(opt.iter)
      
      if(opt.iter=="inf"){ # infinitesimal model
        beta_inf = snp_ldpred2_inf(corr,df_beta,h2=h2_est)
        pred_prs = big_prodVec(snp_fastImputeSimple(genotype),beta_inf,ind.row=ind.test,ind.col=info_snp$"_NUM_ID_")
      }
      if(opt.iter=="grid"){ # grid model : p/h2 needed 
        p_seq = signif(seq_log(1e-4,1,length.out=20),2) # designate length of seq
        h2_seq = round(h2_est*c(0.7,1,1.4),4) # designate length of seq
        grid.param = expand.grid(p=p_seq,h2=h2_seq,sparse=c(FALSE,TRUE))
        grid = grid.param; grid.df = data.frame(PRS=paste0("PRS",1:(20*3*2)),grid); grid.df$PRS = as.character(grid.df$PRS) # designate the length of PRS depending on p&h2
        write.table(grid.df,file=paste0(name.odir,name.test,".LDpred.",opt.iter,".prs.appendix"),row.names=F,quote=F)
        
        beta_grid = snp_ldpred2_grid(corr,df_beta,grid.param,ncores=min(80,ncore))
        pred_prs = big_prodMat(genotype,beta_grid,ind.col=info_snp$"_NUM_ID_")
        colnames(pred_prs) = paste0("PRS",1:(20*3*2)) # designate the length of PRS depending on p&h2
      }
      if(opt.iter=="auto"){ # auto model
        multi_auto = snp_ldpred2_auto(corr,df_beta,h2_init=h2_est,vec_p_init=seq_log(1e-4,0.9,length.out=max(80,ncore)),ncores=min(80,ncore))
        beta_auto = sapply(multi_auto, function(auto){auto$beta_est})
        pred_auto = big_prodMat(genotype,beta_auto,ind.row = ind.test,ind.col = info_snp$`_NUM_ID_`)
        pred_scaled = apply(pred_auto,2,sd) # scale the PRS generated from auto
        beta_auto_rev = rowMeans(beta_auto[,abs(pred_scaled-median(pred_scaled))<3*mad(pred_scaled)])
        pred_prs = big_prodVec(genotype,beta_auto_rev,ind.row=ind.test,ind.col=info_snp$"_NUM_ID_")
      }
      prs.df = data.frame(FID=fam.df$V1,IID=fam.df$V2,PRS=pred_prs)
      write.table(prs.df,file=paste0(name.odir,name.test,".",method,".",opt.iter,".prs"),row.names=F,quote=F)
      
    }
    
  }
  
  if(method=="LASSOsum"){
    
    # Calculate SNP-wise correlations converted from P_val
    ss.df = ss; cor = p2cor(p=ss.df$P,n=ss.df$N,sign=ss.df$BETA) 
    
    # Run Lassosum
    out = lassosum.pipeline(cor=cor,chr=ss.df$CHR,pos=ss.df$BP,A1=ss.df$ALT,A2=ss.df$REF, # A1(ALT)/A2(REF) : A2 is not required but advised
                            ref.bfile=test.bfile,test.bfile=test.bfile,LDblocks="ASN.hg19") # If ref.bfile is unavailable, put test.bfile into ref.bfile.
    v = validate(out) # Validation using the 6th column in .fam to choose lambda and s
    
    # Calculate PRS
    fam.df = fread(paste0(test.bfile,".fam"),data.table=F)
    prs.df = data.frame(FID=fam.df$V1,IID=fam.df$V2,PRS=v$best.pgs)
    write.table(prs.df,file=paste0(name.odir,name.test,".",method,".prs"),row.names=F,quote=F)
    
  }
  
  if(method=="PRScs"){
    
    # Copy the files
    if(!dir.exists(paste0(name.idir,"PRScs-tmp"))){ 
      system(paste("mkdir",paste0(name.idir,"PRScs-tmp")))
      system(paste("mkdir",paste0(name.idir,"PRScs-tmp/PRScs_code")))
      system(paste("mkdir",paste0(name.idir,"PRScs-tmp/ldblk_1kg_eas")))
      system(paste("cp /data4/jyjo/Tutorial/PRS/PRScs-tmp/PRScs_code/*",paste0(name.idir,"PRScs-tmp/PRScs_code/")))
      system(paste("cp /data4/jyjo/Tutorial/PRS/PRScs-tmp/ldblk_1kg_eas/*",paste0(name.idir,"PRScs-tmp/ldblk_1kg_eas/")))
    }
    
    # Rename the colnames of REF/ALT to A2/A1 
    colnames(ss)[colnames(ss)%in%"REF"] = "A2"
    ss.df = ss[,c("SNP","A1","A2","BETA","P")]
    write.table(ss.df,file=paste0(name.idir,"PRScs-tmp/",name.ss,".sumstats.PRScs"),row.names=F,quote=F)
    
    # Calculate posterior SNP effect size estimates for each chromosome
    system(paste("python3",paste0(name.idir,"PRScs-tmp/PRScs_code/PRScs.py"),paste0("--ref_dir=",name.idir,"PRScs-tmp/ldblk_1kg_eas"),
                 paste0("--bim_prefix=",name.idir,name.test),paste0("--sst_file=",name.idir,"PRScs-tmp/",name.ss,".sumstats.PRScs"),
                 paste0("--n_gwas=",max(ss$N)),"--phi=1e-2",paste0("--out_dir=",name.idir,"PRScs-tmp/",name.test)))
    
    # Calculate the PRS using PLINK 
    for(chr.iter in 1:22){ 
      system(paste("plink --bfile",test.bfile,"--score",paste0(name.idir,"PRScs-tmp/",name.test,"_pst_eff_a1_b0.5_phi1e-02_chr",chr.iter,".txt"),"2 4 6 sum",
                   "--out",paste0(name.idir,"PRScs-tmp/",name.test,"_chr",chr.iter),"--threads",min(80,ncore)))
    }
    
    # Combine the chr-separated PRS into one PRS
    tem.df = c()
    for(chr.iter in 1:22){
      tem = fread(paste0(name.idir,"/PRScs-tmp/",name.test,"_chr",chr.iter,".profile"),data.table=F)
      tem.df = cbind(tem.df,tem$SCORESUM)
    }
    prs.vec = rowSums(tem.df)
    prs.df = data.frame(FID=fam.df$V1,IID=fam.df$V2,PRS=prs.vec)
    write.table(prs.df,file=paste0(name.odir,name.test,".",method,".prs"),row.names=F,quote=F)
  }
}


#########################

#### Summary stat 
#### .sumstats should include CHR/BP/SNP/REF/ALT/N/BETA/P columns. 
ss = lapply(1:22, function(i) {
  tmp_dat <- bigreadr::fread2(paste0("/data2/Rex_analysis/BBJ/summary_stat/T2D/T2D.chr", i, ".auto_sex_stratified.txt")) %>%
    mutate(SNP = gsub("_", ":", SNPID), N = 45383+132032) %>%
    dplyr::select(CHR, BP = POS, SNP, REF = Allele1, ALT = Allele2, 
                  BETA = metabeta, N, SE = metase, P = metap)
  return(tmp_dat)
  })
ss <- do.call(rbind, ss)
write.table(ss, "/home2/sjseo/Bohun/SNParray/derived_data/03_04_PRS_BBJ/DM_discovery.sumstats", col.names = T, row.names = F, quote = F)
#ss <- read.table("/home2/sjseo/Bohun/SNParray/derived_data/03_04_PRS_BBJ/DM_discovery.sumstats",
#                 header = T, stringsAsFactors = F)

# summary statistics from discovery dataset : SNUH_BP.sumstats
# test dataset : KARE.bed/.bim/.fam

name.idir = "/home2/sjseo/Bohun/SNParray/derived_data/03_04_PRS_BBJ/"
name.ss = "DM_discovery"
name.test = "02_01_IDmatched"
name.odir = "/home2/sjseo/Bohun/SNParray/derived_data/03_04_PRS_BBJ/"

fam = fread(paste0(paste0(name.idir,name.test),".fam"),data.table=F)
pheno <- read.table("/home2/sjseo/Bohun/SNParray/derived_data/02_01_IDmatch/Matched_ID.txt", header = T)
fam.upd <- fam %>% left_join(., pheno, by = c("V1" = "Methyl")) %>%
  mutate(V6 = ifelse(Pheno == "Normal", 0, 1)) %>%
  select(1:6)
write.table(fam.upd, paste0(paste0(name.idir,name.test),".fam"), col.names = F, row.names = F, quote = F, sep = "\t")

PT.rst = uniPRS.fun(name.idir=name.idir,name.ss=name.ss,name.test=name.test,name.odir=name.odir,method="PT")
CT.rst = uniPRS.fun(name.idir=name.idir,name.ss=name.ss,name.test=name.test,name.odir=name.odir,method="CT")
LDpred.rst = uniPRS.fun(name.idir=name.idir,name.ss=name.ss,name.test=name.test,name.odir=name.odir,method="LDpred")
#LASSOsum.rst = uniPRS.fun(name.idir=name.idir,name.ss=name.ss,name.test=name.test,name.odir=name.odir,method="LASSOsum")
#PRScs.rst = uniPRS.fun(name.idir=name.idir,name.ss=name.ss,name.test=name.test,name.odir=name.odir,method="PRScs")


#### Phenotype
library(dplyr)
library(ggpubr)

pheno <- read.table("/home2/sjseo/Bohun/SNParray/derived_data/02_01_IDmatch/Matched_ID.txt", header = T)
#prs.df <- read.table(paste0(name.odir, "IDmatched.PT.prs"), header = T, stringsAsFactors = F)
#prs.df <- read.table(paste0(name.odir, "IDmatched.CT.prs"), header = T, stringsAsFactors = F)
#prs.df <- read.table(paste0(name.odir, "IDmatched.LDpred.inf.prs"), header = T, stringsAsFactors = F)
prs.df <- read.table(paste0(name.odir, "IDmatched.LASSOsum.prs"), header = T, stringsAsFactors = F)
prs_pheno <- prs.df %>% left_join(., pheno, by = c("IID" = "Methyl"))

#aov_res <- aov(PRS ~ Pheno, data=prs_pheno)
#summary(aov_res)
#posthoc <- TukeyHSD(x=aov_res, 'Pheno', conf.level=0.95)
#posthoc

my_comparisons <- list(c("DMwithAO", "DMwithoutAO"), c("DMwithoutAO", "Normal"), c("DMwithAO", "Normal") )

# size : 400*350
ggplot(data = prs_pheno,  aes(x = Pheno, y = PRS, fill = Pheno)) +
  geom_point() + geom_boxplot() + theme_bw() +
  stat_compare_means(comparisons = my_comparisons)+
  stat_compare_means(method = "anova", label.y = -4) +
  scale_x_discrete(labels = c("AO-exposed\nT2D", "AO-unexposed\nT2D", "Normal")) +
  scale_y_continuous(limits=c(-4.5, 5.8)) +
  scale_fill_manual(values = grey.colors(3, start = 0.3, end = 0.9)) +
  labs(x = NULL) +
  theme(axis.text.x = element_text(color = "black", size =11),
        axis.text.y = element_text(color = "black", size = 11),
        text = element_text(size=13,  family="Arial", color = "black"), 
        legend.position = "none") 