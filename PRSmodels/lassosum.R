#!/usr/bin/env Rscript

###################################################################
# Import packages needed
library(lassosum)
library(data.table)

###############################################################
# parse input arguments
Sys.setlocale("LC_ALL", "C")
options(stringsAsFactors = F)

## Collect arguments
args <- commandArgs(TRUE)

## Default setting when no arguments passed
if(length(args) < 1) {
  args <- c("--help")
}

## Parse arguments (we expect the form --arg=value)
parseArgs <- function(x) strsplit(sub("^--", "", x), "=")
argsDF <- as.data.frame(do.call("rbind", parseArgs(args)))
argsL <- as.list(as.character(argsDF$V2))
names(argsL) <- argsDF$V1

## Check mandatory parameters
if (is.null(argsL$chr)) {
  cat('* Please specify the chromosome --chr\n')
  q(save="no")
} else if (is.null(argsL$bim_file)) {
  cat('* Please specify the directory to the reference PLINK file --bim_file\n')
  q(save="no")
} else if (is.null(argsL$sst_file)) {
  cat('* Please specify the path to the summary statistics with standardized beta using --sst_file\n')
  q(save="no")
} else if (is.null(argsL$LDblocks)) {
  cat('* Please specify the name of LDblocks --LDblocks\n')
  q(save="no")
} else if (is.null(argsL$out_path)) {
  cat('* Please specify the output path\n')
  q(save="no")
} 

print(argsL)

###############################################################
# time calculation
start_time <- Sys.time()

# Create the output file
chr = argsL$chr

# Specify the PLINK file of the reference panel
bfile <- argsL$bim_file
# Read the summary statistics of standardized beta in single variant test
# the standardized beta in single variant test = correlation
ss <- fread(argsL$sst_file)
cor <- ss$Beta

# lassosum only allow -1 < cor < 1
if (sum(abs(cor) >= 1) > 0){
  
  shrink_factor = max(abs(cor)) / 0.9999
  
  cor = cor / shrink_factor
  
}

ss$SNPPos <- sapply(1:length(ss$SNP), function(i) strsplit(ss$SNP[i], "_")[[1]][2])
ss$Chrom <- chr

ref_ld_block <- argsL$LDblocks

print(paste0("Lassosum reference LD blocks: ",ref_ld_block))

if(!ref_ld_block %in% c("EUR.hg19","AFR.hg19","ASN.hg19","EUR.hg38","AFR.hg38","ASN.hg38")){
  stop("LD block does not match blocks defined by Berisa and Pickrell 2015")
}

# train lassosum
out <- lassosum.pipeline(cor=cor,
                         chr=as.numeric(ss$Chrom),
                         pos=as.numeric(ss$SNPPos),
                         A1=ss$A1,
                         A2=ss$A2, # A2 is not required but advised
                         s = c(0.2, 0.5, 0.9, 1),
                         lambda = exp(seq(log(0.0001), log(0.1), length.out = 20)),
                         ref.bfile = bfile, # The reference panel dataset
                         test.bfile = bfile, # We don't have test data here
                         LDblocks = ref_ld_block,
                         exclude.ambiguous = F,
                         destandardize = F,
                         trace = 0)

# perform pseudovalidation
v <- pseudovalidate(out)
lassosum_out <- subset(out, s=v$best.s, lambda=v$best.lambda)

# save estimated beta
sumstats = lassosum_out$sumstats[, c("chr", "pos", "A1", "A2")]
beta = unlist(lassosum_out$beta)
results = data.frame(sumstats, ES = beta) 
results = results[, c("chr", "pos", "A1", "A2", "ES")]

write.table(results,
            argsL$out_path,
            quote = F,
            row.names= F,
            col.names= T,
            sep = "\t",
            append = F)

