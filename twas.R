#!/usr/bin/env Rscript

###################################################################
# parse arguments
###################################################################

Sys.setlocale("LC_ALL", "C")
options(stringsAsFactors = FALSE)

args <- commandArgs(trailingOnly = TRUE)

# Show help if no arguments
if (length(args) < 1) {
  cat("Usage: Rscript twas.R --chr=4 --PRS_methods=PT,lassosum --pheno_file=fileX ... \n")
  q(save = "no")
}

# Parse arguments: --arg=value
parseArgs <- function(x) strsplit(sub("^--", "", x), "=")
argsDF <- as.data.frame(do.call("rbind", parseArgs(args)), stringsAsFactors = FALSE)
argsL <- setNames(as.list(argsDF$V2), argsDF$V1)

# Validate required arguments
required <- c("chrom", "models", "pt", "pheno_file", 
              "grex_dir", "out_dir", "anno_file", "pheno_num")

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
PRS_methods=as.character(argsL$models)
pt=as.character(argsL$pt)
pheno_file=as.character(argsL$pheno_file)
grex_dir=as.character(argsL$grex_dir)
anno_file=as.character(argsL$anno_file)
pheno_num=as.integer(argsL$pheno_num)
out_dir=as.character(argsL$out_dir)

###################################################################
# load dependencies
###################################################################
library(data.table)
library(ACAT)

######################################################################
# load helper functions
######################################################################

do_twas_gene <- function(gene,model,pheno_name,verbose=F){
  if(verbose){
    print(c(gene,model,pheno_name))
  }
  
  tmp <- dat_chr[dat_chr$model==model & dat_chr$gene==gene,]
  if(nrow(tmp)>0){
    tmp <- merge(x=tmp,y=pheno,by.x="id",by.y="SampID",all.x=T)
    
    tmp$grex_comb_lai <- tmp$aspPS_anc0+tmp$aspPS_anc1
    tmp$pheno <- tmp[,which(names(tmp)==pheno_name)]

    # determine if phenotype is binary otherwise fit linear regression
    unique_vals <- unique(tmp$pheno)
    unique_vals <- unique_vals[order(unique_vals)]
    if(length(unique_vals) == 2 && isTRUE(all.equal(sort(unique_vals), c(0, 1)))){

      # fit test linear model to extract column names
      y <- sample(0:1,100,replace=T);x <- rnorm(100)
      fit <- glm(y~x,family = "binomial")
      tmp_names <- colnames(t(summary(fit)$coef[2,]))

      # initialize results
      res <- {}
      
      # Fit a logistic model
      fit <- glm(pheno~grex_comb_lai,dat=tmp,family="binomial")
      if(dim(summary(fit)$coef)[1]>1){
        res <- rbind(res,as.data.frame(t(summary(fit)$coef[2,])))
      }else{
        tmp2 <- data.frame(matrix(nrow=1,ncol=4))
        names(tmp2) <- tmp_names
        res <- rbind(res,tmp2)
      }
      
      fit <- glm(pheno~aspPS_anc0,dat=tmp,family="binomial")
      if(dim(summary(fit)$coef)[1]>1){
        res <- rbind(res,as.data.frame(t(summary(fit)$coef[2,])))
      }else{
        tmp2 <- data.frame(matrix(nrow=1,ncol=4))
        names(tmp2) <- tmp_names
        res <- rbind(res,tmp2)
      }
      
      fit <- glm(pheno~aspPS_anc1,dat=tmp,family="binomial")
      if(dim(summary(fit)$coef)[1]>1){
        res <- rbind(res,as.data.frame(t(summary(fit)$coef[2,])))
      }else{
        tmp2 <- data.frame(matrix(nrow=1,ncol=4))
        names(tmp2) <- tmp_names
        res <- rbind(res,tmp2)
      }
      
      fit <- glm(pheno~PS_anc0,dat=tmp,family="binomial")
      if(dim(summary(fit)$coef)[1]>1){
        res <- rbind(res,as.data.frame(t(summary(fit)$coef[2,])))
      }else{
        tmp2 <- data.frame(matrix(nrow=1,ncol=4))
        names(tmp2) <- tmp_names
        res <- rbind(res,tmp2)
      }
      
      fit <- glm(pheno~PS_anc1,dat=tmp,family="binomial")
      if(dim(summary(fit)$coef)[1]>1){
        res <- rbind(res,as.data.frame(t(summary(fit)$coef[2,])))
      }else{
        tmp2 <- data.frame(matrix(nrow=1,ncol=4))
        names(tmp2) <- tmp_names
        res <- rbind(res,tmp2)
      }
      
      res$method <- c("casPS",
                      "aspPS-anc0",
                      "aspPS-anc1",
                      "PS-anc0",
                      "PS-anc1")
      res$regression <- "logit"

    }else{

      # fit test linear model to extract column names
      y <- rnorm(100);x <- rnorm(100)
      fit <- lm(y~x)
      tmp_names <- colnames(t(summary(fit)$coef[2,]))

      # initialize results
      res <- {}
      
      # Fit a linear model
      fit <- lm(pheno~grex_comb_lai,dat=tmp)
      if(dim(summary(fit)$coef)[1]>1){
        res <- rbind(res,as.data.frame(t(summary(fit)$coef[2,])))
      }else{
        tmp2 <- data.frame(matrix(nrow=1,ncol=4))
        names(tmp2) <- tmp_names
        res <- rbind(res,tmp2)
      }
      
      fit <- lm(pheno~aspPS_anc0,dat=tmp)
      if(dim(summary(fit)$coef)[1]>1){
        res <- rbind(res,as.data.frame(t(summary(fit)$coef[2,])))
      }else{
        tmp2 <- data.frame(matrix(nrow=1,ncol=4))
        names(tmp2) <- tmp_names
        res <- rbind(res,tmp2)
      }
      
      fit <- lm(pheno~aspPS_anc1,dat=tmp)
      if(dim(summary(fit)$coef)[1]>1){
        res <- rbind(res,as.data.frame(t(summary(fit)$coef[2,])))
      }else{
        tmp2 <- data.frame(matrix(nrow=1,ncol=4))
        names(tmp2) <- tmp_names
        res <- rbind(res,tmp2)
      }
      
      fit <- lm(pheno~PS_anc0,dat=tmp)
      if(dim(summary(fit)$coef)[1]>1){
        res <- rbind(res,as.data.frame(t(summary(fit)$coef[2,])))
      }else{
        tmp2 <- data.frame(matrix(nrow=1,ncol=4))
        names(tmp2) <- tmp_names
        res <- rbind(res,tmp2)
      }
      
      fit <- lm(pheno~PS_anc1,dat=tmp)
      if(dim(summary(fit)$coef)[1]>1){
        res <- rbind(res,as.data.frame(t(summary(fit)$coef[2,])))
      }else{
        tmp2 <- data.frame(matrix(nrow=1,ncol=4))
        names(tmp2) <- tmp_names
        res <- rbind(res,tmp2)
      }
      
      res$method <- c("casPS",
                      "aspPS-anc0",
                      "aspPS-anc1",
                      "PS-anc0",
                      "PS-anc1")
      res$regression <- "linear"
    }
    
    res$gene=gene
    res$model=model
    res$pheno=pheno_name
    return(res)
  }else{
    cat(paste0("No rows in GReX file, skipping: ",model," for ", gene,"\n"))
  }  
}

######################################################################
# analyze data and save results:
######################################################################

# Check if the output directory exists
if (!dir.exists(out_dir)) {
  # If it doesn't exist, create it
  dir.create(out_dir)
}

# Read in phenotype and gene annotation data
pheno <- data.frame(fread(pheno_file))
anno <- data.frame(fread(anno_file))

# Subset anno file to chromosome of interest
anno <- anno[anno$CHROM==chrom,]
genes <- anno$TargetID

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
fnames_methods <- paste0("chr",chrom,"_",fnames_methods)

# Check if all required files exist in imputed grex directory
file_list <- list.files(grex_dir)
check <- which(fnames_methods %in% file_list)
if (length(check) != length(fnames_methods)) {
  missing <- fnames_methods[!fnames_methods %in% file_list]
  stop(paste0("GReX files not found for all PRS methods: ", paste(missing, collapse = ", ")))
}

# Read in imputed grex files across all PRS models
dat <- data.frame()
for (f in fnames_methods) {
  file_path <- file.path(grex_dir, f)
  tmp <- fread(file_path)
  if (nrow(tmp) > 0) {
    dat <- rbind(dat, tmp)
  } else {
    stop(paste0("Empty GReX file: ", f))
  }
}
dat_chr <- data.frame(dat)

dat_chr$chr <- chrom
dat_chr <- dat_chr[dat_chr$gene %in% genes,]

if(nrow(dat_chr)>0){
  methods <- unique(dat_chr$model)
  
  pheno_vec <- names(pheno)[-1]
  pheno_name <- pheno_vec[pheno_num]
  genes <- unique(dat_chr$gene)
  print(paste0("Total number of genes available for TWAS on this chromosome: ",length(genes)))

  start.time <- Sys.time()
  res <- lapply(unique(dat_chr$gene),FUN=function(x){
    out <- {}
    for(method in methods){
        tmp <- do_twas_gene(gene=x,model=method,pheno_name=pheno_name,verbose = F)
        if(!is.null(tmp)){
          out <- rbind(out,tmp)
        }
    }
    return(out)
  })
  
  end.time <- Sys.time()
  run.time <- difftime(end.time,start.time,units="mins")
  cat("TWAS regression models completed.\n")
  print("Total run time: ")
  print(run.time)
  
  dat <- do.call(rbind, res)
  dat$chr <- chrom
  
}else{
  print("No genes with imputed GRex found for genes in annotation file for this chromosome; terminating process")
}

cat("Completed TWAS p-value calculation. Moving onto to stage I and II pvalue combination with ACAT. \n")


#### ACAT p-value combination

calc_p <- function(gene,pheno_name){
  tmp <- dat[dat$gene==gene & dat$pheno==pheno_name,]

  if (tmp$regression[1]=="linear"){
      tmp <- tmp[!is.na(tmp$`Pr(>|t|)`),]
      if(nrow(tmp)>3){
      out <- {}
      out <- rbind(out,c("casPS", ACAT(tmp$`Pr(>|t|)`[tmp$method=="casPS"])))
      out <- rbind(out,c("aspPS", ACAT(tmp$`Pr(>|t|)`[tmp$method%in%c("aspPS-anc0","aspPS-anc1")])))
      out <- rbind(out,c("PS-anc0", ACAT(tmp$`Pr(>|t|)`[tmp$method%in%c("PS-anc0")])))
      out <- rbind(out,c("PS-anc1", ACAT(tmp$`Pr(>|t|)`[tmp$method%in%c("PS-anc1")])))
      out <- rbind(out,c("CADET", ACAT(tmp$`Pr(>|t|)`[tmp$method%in%c("casPS","aspPS-anc0","aspPS-anc1","PS-anc0","PS-anc1")])))
      out <- data.frame(out)
      out$gene <- gene
      out$pheno <- pheno_name
      names(out)[1:2] <- c("ACAT Approach","p")
    }
  }else if(tmp$regression[1]=="logit"){
      tmp <- tmp[!is.na(tmp$`Pr(>|z|)`),]
      if(nrow(tmp)>3){
      out <- {}
      out <- rbind(out,c("casPS", ACAT(tmp$`Pr(>|z|)`[tmp$method=="casPS"])))
      out <- rbind(out,c("aspPS", ACAT(tmp$`Pr(>|z|)`[tmp$method%in%c("aspPS-anc0","aspPS-anc1")])))
      out <- rbind(out,c("PS-anc0", ACAT(tmp$`Pr(>|z|)`[tmp$method%in%c("PS-anc0")])))
      out <- rbind(out,c("PS-anc1", ACAT(tmp$`Pr(>|z|)`[tmp$method%in%c("PS-anc1")])))
      out <- rbind(out,c("CADET", ACAT(tmp$`Pr(>|z|)`[tmp$method%in%c("casPS","aspPS-anc0","aspPS-anc1","PS-anc0","PS-anc1")])))
      out <- data.frame(out)
      out$gene <- gene
      out$pheno <- pheno_name
      names(out)[1:2] <- c("ACAT Approach","p")
      
      }
  }
  return(out)
}

# aggregate final results
gene_list <- unique(dat$gene)
out <- lapply(gene_list,FUN=calc_p,pheno_name=pheno_name)
res <- do.call(rbind, out)

write.table(res,file=paste0(out_dir,"/chr",chrom,"_twas_",pheno_name,".txt"),
            sep="\t",col.names = T,row.names = F,quote=F)

cat("Completed TWAS p-value combination. Analysis complete. \n")

# print memory usage
gc()










