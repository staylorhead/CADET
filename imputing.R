
#!/usr/bin/env Rscript

##############################################################################
# load dependencies
##############################################################################
library(data.table)

##############################################################################
# parse arguments
##############################################################################

Sys.setlocale("LC_ALL", "C")
options(stringsAsFactors = FALSE)

args <- commandArgs(trailingOnly = TRUE)

# Show help if no arguments
if (length(args) < 1) {
  cat("Usage: Rscript imputing.R --chr=4 --anc0_models_dir=dir0 --anc1_models_dir=dir1 ... \n")
  q(save = "no")
}

# Parse arguments: --arg=value
parseArgs <- function(x) strsplit(sub("^--", "", x), "=")
argsDF <- as.data.frame(do.call("rbind", parseArgs(args)), stringsAsFactors = FALSE)
argsL <- setNames(as.list(argsDF$V2), argsDF$V1)

# Validate required arguments
required <- c("chrom", "anc0_models_dir", "anc1_models_dir", "vcf_file", 
              "anno_file", "models", "pt", "maf_anc0", "maf_anc1", "out_dir")

missing <- setdiff(required, names(argsL))
if (length(missing) > 0) {
  cat("* Missing required arguments:\n", paste0("--", missing, collapse = "\n"), "\n")
  q(save = "no")
}

# Print parsed arguments 
cat("Parsed arguments:\n")
print(argsL)

# convert types
chrom=as.integer(argsL$chrom)
anc0_models_dir=as.character(argsL$anc0_models_dir)
anc1_models_dir=as.character(argsL$anc1_models_dir)
PRS_methods=as.character(argsL$models)
pt=as.character(argsL$pt)
vcf_file=as.character(argsL$vcf_file)
anno_file=as.character(argsL$anno_file)
maf_anc0=as.character(argsL$maf_anc0)
maf_anc1=as.character(argsL$maf_anc1)
out_dir=as.character(argsL$out_dir)

##############################################################################
# load helper functions
##############################################################################

# function to extract hap data for admixed samples
extract_haplotype_adm <- function(vcf, hapnum, data){
  if(data=="genotype"){
    hap <- apply(vcf,2,function(x){
      tmp <- strsplit(x,":")
      geno <- unlist(lapply(tmp, `[[`, 1))
      hap <- as.numeric(substr(geno,2*hapnum-1,2*hapnum-1))
      return(hap)
    })
  } else if(data=="ancestry"){
    hap <- apply(vcf,2,function(x){
      tmp <- strsplit(x,":")
      if(hapnum==1){
        hap <- unlist(lapply(tmp, `[[`, 2))
      }else if (hapnum==2){
        hap <- unlist(lapply(tmp, `[[`, 3))
      }
      return(hap)
    })
  }
  return(hap)
}

# function to standardize hap data
std_test_hap <- function(i,hap_dat,ancestry_dat,snp_info){
  hap <- as.vector(hap_dat[,i])
  ancestry <- as.vector(ancestry_dat[,i])
  anc0_idx <- which(ancestry=="1")
  anc1_idx <- which(ancestry=="0")

  maf_anc0 <- snp_info$maf_anc0 
  maf_anc1 <- snp_info$maf_anc1
  
  if(length(anc0_idx)>1){
    tmp <- sapply(anc0_idx,function(x){
      (hap[x]-2*maf_anc0[x])/sqrt(2*maf_anc0[x]*(1-maf_anc0[x]))
    })
    hap[anc0_idx] <- tmp
  }
  if(length(anc1_idx)>1){
    tmp <- sapply(anc1_idx,function(x){
      (hap[x]-2*maf_anc1[x])/sqrt(2*maf_anc1[x]*(1-maf_anc1[x]))
    })
    hap[anc1_idx] <- tmp
  }
  return(hap)
}

# function to impute grex using LAI
impute_grex_anc <- function(wt,hap1_geno_test,hap2_geno_test,
                            hap1_anc_test,hap2_anc_test,n_samp_test,reference,
                            anc_std=F){
  
    if(nrow(wt)>0){
      wt <- wt[!wt$ES==0,]
      if(nrow(wt)>0){
        snpID_wt <- wt$SNP
        
        if(nrow(wt)>1){
          X1 <- hap1_geno_test[snpID_wt,]
          X2 <- hap2_geno_test[snpID_wt,]
          
          causal_anc0 <- hap1_anc_test[snpID_wt,]
          causal_anc1 <- hap2_anc_test[snpID_wt,]
        }else{
          X1 <- matrix(hap1_geno_test[snpID_wt,],nrow=1)
          X2 <- matrix(hap2_geno_test[snpID_wt,],nrow=1)
          
          causal_anc0 <- matrix(hap1_anc_test[snpID_wt,],nrow=1)
          causal_anc1 <- matrix(hap2_anc_test[snpID_wt,],nrow=1)
        }
        
        Gamma1 <- matrix(0,nrow=nrow(X1),ncol=ncol(X1))
        Gamma2 <- matrix(0,nrow=nrow(X1),ncol=ncol(X1))
        
        Gamma1[causal_anc0==reference] <- 1
        Gamma2[causal_anc1==reference] <- 1
        
        # some values will be 0 because they may be, in a indiv, of another
        # ancestry and that ancestry doesn't have MAF for that snp
        if(anc_std){
          X1[is.na(X1)] <- 0
          X2[is.na(X2)] <- 0
        }
        
        beta <- wt$ES
        grex <- t(X1*Gamma1+X2*Gamma2)%*%beta
        
      }
      
    }else{
      grex <- matrix(NA,nrow=n_samp_test,ncol=1)
    }

  
  return(grex=grex)
}

# function to impute grex using standard approach
impute_grex <- function(wt,geno_dip_test,n_samp_test){
    if(nrow(wt)>0){
      wt <- wt[!wt$ES==0,]
      if(nrow(wt)>0){
        snpID_wt <- wt$SNP
        if(nrow(wt)>1){
          X <- geno_dip_test[snpID_wt,]
          beta <- wt$ES
          grex <- t(X)%*%beta
        }else{
          X <- geno_dip_test[snpID_wt,]
          beta <- wt$ES
          grex <- X*beta
        }
        
      }
      
    }else{
      grex <- matrix(NA,nrow=n_samp_test,ncol=1)
    }
  return(grex)
}

# function to impute grex for a single gene
impute_grex_gene <- function(gene, vcf, anc0_models, anc1_models){
  print(gene)
  
  weight_anc0 <- anc0_models[anc0_models$TargetID==gene,]
  weight_anc1 <- anc1_models[anc1_models$TargetID==gene,]
  
  weight_anc0$SNP <- paste0(weight_anc0$CHROM,"_",weight_anc0$POS,"_",weight_anc0$A2,"_",weight_anc0$A1)
  weight_anc1$SNP <- paste0(weight_anc1$CHROM,"_",weight_anc1$POS,"_",weight_anc1$A2,"_",weight_anc1$A1)
  weight_anc0$SNP2 <- paste0(weight_anc0$CHROM,"_",weight_anc0$POS,"_",weight_anc0$A1,"_",weight_anc0$A2)
  weight_anc1$SNP2 <- paste0(weight_anc1$CHROM,"_",weight_anc1$POS,"_",weight_anc1$A1,"_",weight_anc1$A2)

  # Subset and flip ES if SNP is in reverse orientation (SNP2)
  weight_anc0$keep <- weight_anc0$SNP %in% vcf$ID
  weight_anc0$flip <- weight_anc0$SNP2 %in% vcf$ID
  weight_anc1$keep <- weight_anc1$SNP %in% vcf$ID
  weight_anc1$flip <- weight_anc1$SNP2 %in% vcf$ID

  # Keep only those present in either orientation
  weight_anc0 <- weight_anc0[weight_anc0$keep | weight_anc0$flip, ]
  weight_anc1 <- weight_anc1[weight_anc1$keep | weight_anc1$flip, ]

  # Flip ES where SNP2 matches VCF
  weight_anc0$ES[weight_anc0$flip] <- -weight_anc0$ES[weight_anc0$flip]
  weight_anc1$ES[weight_anc1$flip] <- -weight_anc1$ES[weight_anc1$flip]

  # Set the SNP ID to the correct matching form (SNP or SNP2)
  weight_anc0$SNP <- ifelse(weight_anc0$flip, weight_anc0$SNP2, weight_anc0$SNP)
  weight_anc1$SNP <- ifelse(weight_anc1$flip, weight_anc1$SNP2, weight_anc1$SNP)

  # Clean up and sort by position
  weight_anc0 <- weight_anc0[ , !(names(weight_anc0) %in% c("SNP2", "keep", "flip"))]
  weight_anc0 <- weight_anc0[order(weight_anc0$POS), ]
  weight_anc1 <- weight_anc1[ , !(names(weight_anc1) %in% c("SNP2", "keep", "flip"))]
  weight_anc1 <- weight_anc1[order(weight_anc1$POS), ]

  # remove SNPs for each ancestry that are missing from ancestry-specific MAF dadta
  weight_anc0 <- weight_anc0[weight_anc0$SNP %in% ss_data_anc0$ID,]
  weight_anc1 <- weight_anc1[weight_anc1$SNP %in% ss_data_anc1$ID,]

  snp_comb <- unique(c(weight_anc0$SNP,weight_anc1$SNP))
  #sum(weight_anc0$SNP %in% weight_anc1$SNP)
  
  vcf_comb <- vcf[vcf$ID %in% snp_comb,]
  snp_info <- vcf_comb[,1:9]
  rownames(snp_info) <- snp_info$ID
  orig_order <- rownames(snp_info)
    
  weight_anc0 <- merge(x=weight_anc0,y=ss_data_anc0,
                      by.x="SNP",by.y="ID",all.x=T,all.y=F)
  
  weight_anc1 <- merge(x=weight_anc1,y=ss_data_anc1,
                      by.x="SNP",by.y="ID",all.x=T,all.y=F)
  
  snp_info_2 <- merge(x=snp_info,y=weight_anc0[,c("SNP","MAF")],by.x="ID",by.y="SNP",all.x=T,all.y=F)
  names(snp_info_2)[ncol(snp_info_2)] <- "maf_anc0"
  
  snp_info <- merge(x=snp_info_2,y=weight_anc1[,c("SNP","MAF")],by.x="ID",by.y="SNP",all.x=T,all.y=F)
  names(snp_info)[ncol(snp_info)] <- "maf_anc1"

  rownames(snp_info) <- snp_info$ID
  snp_info <- snp_info[orig_order,]
  
  vcf_comb <- as.matrix(vcf_comb[,-1:-9])
  
  comb_ancestry_1 <- extract_haplotype_adm(vcf_comb, hapnum=1, data="ancestry")
  comb_ancestry_2 <- extract_haplotype_adm(vcf_comb, hapnum=2, data="ancestry")
  
  comb_hap_1 <- extract_haplotype_adm(vcf_comb, hapnum=1, data="genotype")
  comb_hap_2 <- extract_haplotype_adm(vcf_comb, hapnum=2, data="genotype")
  
  if(is.null(dim(comb_hap_1))){
    comb_hap_1 <- t(as.data.frame(comb_hap_1))
  }
  if(is.null(dim(comb_hap_2))){
    comb_hap_2 <- t(as.data.frame(comb_hap_2))
  }
  if(is.null(dim(comb_ancestry_1))){
    comb_ancestry_1 <- t(as.data.frame(comb_ancestry_1))
  }
  if(is.null(dim(comb_ancestry_2))){
    comb_ancestry_2 <- t(as.data.frame(comb_ancestry_2))
  }
  
  row.names(comb_hap_1) <- snp_info$ID
  row.names(comb_hap_2) <- snp_info$ID
  row.names(comb_ancestry_1) <- snp_info$ID
  row.names(comb_ancestry_2) <- snp_info$ID
  
  n_samp_test <- ncol(vcf_comb)
  
  std_hap1 <- sapply(1:n_samp_test,FUN=std_test_hap,
                     hap_dat=comb_hap_1,
                     ancestry_dat=comb_ancestry_1,
                     snp_info=snp_info)
  row.names(std_hap1) <- snp_info$ID
  std_hap2 <- sapply(1:n_samp_test,FUN=std_test_hap,
                     hap_dat=comb_hap_2,
                     ancestry_dat=comb_ancestry_2,
                     snp_info=snp_info)
  row.names(std_hap2) <- snp_info$ID
  
  grex_anc0 <- impute_grex_anc(wt=weight_anc0,
                              hap1_geno_test=comb_hap_1,
                              hap2_geno_test=comb_hap_2,
                              hap1_anc_test=comb_ancestry_1,
                              hap2_anc_test=comb_ancestry_2,
                              n_samp_test=n_samp_test,
                              reference="0",
                              anc_std = F)
  grex_anc1 <- impute_grex_anc(wt=weight_anc1,
                              hap1_geno_test=comb_hap_1,
                              hap2_geno_test=comb_hap_2,
                              hap1_anc_test=comb_ancestry_1,
                              hap2_anc_test=comb_ancestry_2,
                              n_samp_test=n_samp_test,
                              reference="1",
                              anc_std = F)
  comb_geno <- comb_hap_1+comb_hap_2
  
  grex_std_anc0_ss <- impute_grex(wt=weight_anc0,
                                 geno_dip_test=comb_geno,
                                 n_samp_test=n_samp_test)
  grex_std_anc1_ss <- impute_grex(wt=weight_anc1,
                                 geno_dip_test=comb_geno,
                                 n_samp_test=n_samp_test)
  
  out <- data.frame(gene=gene,
                    grex_anc0=grex_anc0,
                    grex_anc1=grex_anc1,
                    grex_std_anc0=grex_std_anc0_ss,
                    grex_std_anc1=grex_std_anc1_ss,
                    id=colnames(vcf_comb))
  return(out)
}

# function to impute grex for all genes for a single PRS model
impute_grex_method <- function(model, vcf, chrom){
    
  anc0_models <- data.frame(fread(paste0(anc0_models_dir,"/",model),fill=T))
  anc1_models <- data.frame(fread(paste0(anc1_models_dir,"/",model),fill=T))

  # then call and sapply over impute_grex_gene  
  genes_anc0 <- unique(anc0_models$TargetID)
  genes_anc1 <- unique(anc1_models$TargetID)
  olap_genes <- intersect(genes_anc0, genes_anc1)

  print(paste0("Beginning imputation for PRS model: ", model))
  print(paste0("Number of genes with trained GReX models in ancestry 1: ",length(genes_anc0)))
  print(paste0("Number of genes with trained GReX models in ancestry 2: ",length(genes_anc1)))
  print(paste0("Number of genes with trained GReX models in both ancestries: ",length(olap_genes)))
  
  res_list <- lapply(olap_genes, impute_grex_gene, vcf=vcf, anc0_models=anc0_models, anc1_models=anc1_models)
  res <- do.call(rbind, res_list)
  res$model <- model

  # rename columns to be more interpretable
  colnames(res)[colnames(res)=="grex_anc0"] <- "aspPS_anc0"
  colnames(res)[colnames(res)=="grex_anc1"] <- "aspPS_anc1"
  colnames(res)[colnames(res)=="grex_std_anc0"] <- "PS_anc0"
  colnames(res)[colnames(res)=="grex_std_anc1"] <- "PS_anc1"

  write.table(res,paste0(out_dir,"/chr",chrom,"_",model),sep="\t",quote=F,row.names=F,col.names=T)
  
  return(res)
}

##############################################################################
# begin  code to impute grex
##############################################################################

# Check if the output directory exists
if (!dir.exists(out_dir)) {
  # If it doesn't exist, create it
  dir.create(out_dir)
}

# Expand PRS methods
methods_unlist <- strsplit(PRS_methods,",")[[1]] 

# Generate expected GReX model filenames
fnames_methods <- {}
if("PT" %in% methods_unlist){
  pt_vals <- strsplit(pt,",")[[1]] 
  fnames_methods <- c(fnames_methods,paste0("P",pt_vals,".txt"))
}
if("lassosum" %in% methods_unlist){
  fnames_methods <- c(fnames_methods, "lassosum.txt")
}

# Check if all required files exist in anc0 model directory
file_list_anc0 <- list.files(anc0_models_dir)
check_anc0 <- which(fnames_methods %in% file_list_anc0)
if (length(check_anc0) != length(fnames_methods)) {
  missing <- fnames_methods[!fnames_methods %in% file_list_anc0]
  stop(paste0("GReX model(s) not found for all PRS methods for ancestry 1: ", paste(missing, collapse = ", ")))
}

# Check if all required files exist in anc1 model directory
file_list_anc1 <- list.files(anc1_models_dir)
check_anc1 <- which(fnames_methods %in% file_list_anc1)
if (length(check_anc1) != length(fnames_methods)) {
  missing <- fnames_methods[!fnames_methods %in% file_list_anc1]
  stop(paste0("GReX model(s) not found for all PRS methods for ancestry 2: ", paste(missing, collapse = ", ")))
}

# Read in gene annotation file and vcf file with LA information
anno <- data.frame(fread(anno_file))

# Subset to chromosome of interest
anno <- anno[anno$CHROM==chrom,]

# read in vcf file
vcf <- data.frame(fread(vcf_file))
colnames(vcf)[1] <- "#CHROM"
# subset to chomosome
vcf <- vcf[vcf$`#CHROM`==chrom,]

vcf$ID <- paste0(vcf$`#CHROM`,"_",vcf$POS,"_",vcf$REF,"_",vcf$ALT)
vcf$INFO <- paste0(vcf$`#CHROM`,"_",vcf$POS,"_",vcf$ALT,"_",vcf$REF)

# read in MAF data
ss_data_anc0 <- fread(maf_anc0)
ss_data_anc1 <- fread(maf_anc1)

# take 1-MAF if alleles are swapped in vcf
split_id <- strsplit(ss_data_anc0$ID, "_")
ID_flip <- sapply(split_id, function(x) paste(x[1], x[2], x[4], x[3], sep = "_"))
flip_idx <- ID_flip %in% vcf$ID
ss_data_anc0$MAF[flip_idx] <- 1 - ss_data_anc0$MAF[flip_idx]
ss_data_anc0$ID[flip_idx] <- ID_flip[flip_idx]

split_id <- strsplit(ss_data_anc1$ID, "_")
ID_flip <- sapply(split_id, function(x) paste(x[1], x[2], x[4], x[3], sep = "_"))
flip_idx <- ID_flip %in% vcf$ID
ss_data_anc1$MAF[flip_idx] <- 1 - ss_data_anc1$MAF[flip_idx]
ss_data_anc1$ID[flip_idx] <- ID_flip[flip_idx]

# begin imputation
start.time <- Sys.time()
out <- lapply(fnames_methods,impute_grex_method, vcf=vcf, chrom=chrom)
end.time <- Sys.time()
time_cadet <- difftime(end.time,start.time,units="mins")

cat("Completed GReX imputation.\n")

print("Total computation time:")
print(time_cadet)
print(gc())

