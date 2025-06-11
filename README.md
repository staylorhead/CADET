# CADET: Enhanced transcriptome-wide association analyses in admixed samples using eQTL summary data 

## Introduction

<p align="center">
  <img src="images/workflow.png" alt="Figure description" width="600"/>
</p>

CADET is a tool for conducting powerful transcriptome-wide association studies (TWAS) in admixed populations by integrating local ancestry (LA) information and summary-level eQTL data from ancestrally diverse reference panels. CADET predicts local ancestry-aware genetically regulated gene expression (GReX) by combining polygenic risk score (PRS) models with eQTL summary statistics, without requiring individual-level genotype or expression data. This enables downstream TWAS to identify gene-phenotype associations while accounting for possible ancestry-specific genetic regulation.

CADET leverages local ancestry to capture ancestry-dependent architecture of gene expression, improving the relevance and accuracy of TWAS in admixed populations. With this software, you can perform the following core steps:
1.	Train PRS models of gene expression using eQTL summary statistics via:
    - Pruning and thresholding (P+T), or
    - Lassosum
2.	Impute GReX (both LA-aware and LA-unaware) for target admixed individuals
3.	Test gene-phenotype associations via TWAS using the imputed GReX

If you use CADET in your research, please cite the following article:

> Head, S. T., Dai, Q., Schildkraut, J., Cutler, D. J., Yang, J., & Epstein, M. P. (2024). *CADET: Enhanced transcriptome-wide association analyses in admixed samples using eQTL summary data.* bioRxiv. https://doi.org/10.1101/2024.10.21.619441

## Installation 

First, download all software and example files using:

```bash
git clone https://github.com/staylorhead/CADET.git
```

### Dependencies
CADET depends opon successful prior installation of:
- [R](https://www.r-project.org/) and the following packages
    - [lassosum](https://github.com/tshmak/lassosum)
    - [fdrtool](https://cran.r-project.org/web/packages/fdrtool/index.html)
    - [data.table](https://cran.r-project.org/web/packages/data.table/index.html)
    - [ACAT](https://github.com/yaowuliu/ACAT)
- [BGZIP](http://www.htslib.org/doc/bgzip.html)
- [TABIX](http://www.htslib.org/doc/tabix.html)
- [PLINK 1.9](https://www.cog-genomics.org/plink/)
- Python modules/libraries
    - [pandas 1.4.4](https://pandas.pydata.org)
    - [scipy 1.7.3](https://scipy.org)
    - [numpy 1.21.5](https://numpy.org)
    - [pysam 0.19.1](https://pysam.readthedocs.io/en/latest/api.html) 

Here is some example code to create an Python environment for CADET:

```bash
# create environment 
conda create --name cadet python=3.9 pandas=1.4.4 numpy=1.21.5 scipy=1.7.3 pip
# deactivate conda environment
conda deactivate

 # activate the environment
conda activate cadet

# install pysam
pip install pysam

# deactivate the environment
conda deactivate
```

### Input Files

### Example Code

```bash
module load R
module load tabix
module load plink # v1.9
conda activate cadet

# set the location of the CADET source directory
DIR="/home/sfisch3/CADET"

# set the location of user input files
VCF="/home/sfisch3/CADET/Example/flare_geno_anc.vcf"
PHENO="/home/sfisch3/CADET/Example/twas_phenotypes.txt"
ANNO="/home/sfisch3/CADET/Example/anno.txt"
EQTL_SS_ANC0="${DIR}/Example/eQTLSumStatsAnc0.txt"
EQTL_SS_ANC1="${DIR}/Example/eQTLSumStatsAnc1.txt"
MAF_ANC0=${DIR}/Example/mafAnc0.txt
MAF_ANC1=${DIR}/Example/mafAnc1.txt

# first create binary files of target vcf
plink \
  --vcf ${VCF} \
  --double-id \
  --make-bed \
  --keep-allele-order \
  --out ${DIR}/Example/geno

# train eQTL weights in ancestry 0, e.g., AFR
python3 ${DIR}/training.py \
--anno_file=${ANNO} \
--geno_dir=${DIR}/Example/geno \
--out_dir=${DIR}/Example/Output/Anc0_grex_models \
--sst_file=${EQTL_SS_ANC0} \
--lassosum_LD_block="AFR.hg38" \
--r2=0.99 \
--window=1000000 \
--models=PT,lassosum \
--pt=0.001,0.05 \
--thread=1 \
--script_dir=${DIR} \
--seed=123

# train eQTL weights in ancestry 1, e.g., EUR
python3 ${DIR}/training.py \
--anno_file=${ANNO} \
--geno_dir=${DIR}/Example/geno \
--out_dir=${DIR}/Example/Output/Anc1_grex_models \
--sst_file=${EQTL_SS_ANC1} \
--lassosum_LD_block="EUR.hg38" \
--r2=0.99 \
--window=1000000 \
--models=PT,lassosum \
--pt=0.001,0.05 \
--thread=1 \
--script_dir=${DIR} \
--seed=123

# impute GReX for genes on a single chromosome
# can loop over all 22 autosomes
Rscript ${DIR}/imputing.R \
--chrom=4 \
--anc0_models_dir=${DIR}/Example/Output/Anc0_grex_models \
--anc1_models_dir=${DIR}/Example/Output/Anc1_grex_models \
--models=PT,lassosum \
--pt=0.001,0.05 \
--vcf_file=${VCF} \
--anno_file=${ANNO} \
--maf_anc0=${MAF_ANC0} \
--maf_anc1=${MAF_ANC1} \
--out_dir=${DIR}/Example/Output/Imputed_grex

# perform TWAS
# can loop over all 22 autosomes and all phenotypes, e.g.:
for pheno_num in 1 2; do
  Rscript ${DIR}/twas.R \
    --chrom=4 \
    --pheno_num=${pheno_num} \
    --models=PT,lassosum \
    --pt=0.001,0.05 \
    --pheno_file=${PHENO} \
    --grex_dir=${DIR}/Example/Output/Imputed_grex \
    --anno_file=${ANNO} \
    --out_dir=${DIR}/Example/Output/TWAS
done
```
