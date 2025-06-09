# CADET: Enhanced transcriptome-wide association analyses in admixed samples using eQTL summary data 

## Introduction

CADET is a tool that enables powerful transcriptome-wide association study (TWAS) of admixed cohorts leveraging the local-ancestry (LA) information of the cohort along with summary-level eQTL data from reference panels of different ancestral groups. CADET combines multiple polygenic risk score models with the summary-level eQTL reference data to predict LA-aware genetically-regulated gene expression (GReX) in target admixed samples. This gene expression can then be used in TWAS to test for gene-level association with a given phenotype. This software circumvents the need for indiviudal-level genotype and gene expression data and leverage local ancestry information to allow for LA-dependent genetic architecture of gene expression.

There are primary steps performed in CADET:
1. Training PRS models of gene expression using eQTL summary statistics using
    - Pruning + thresholding (P+T), or
    - Lassosum
2. Imputing LA-aware and LA-unaware vectors of gene expression in target admixed subjects
3. Testing for TWAS association with trait of interest

For questions or issues related to this software, please contact Taylor Head (<sthead@mdanderson.org>).

## Installation 

First, download all software and example files using:

```bash
git clone https://github.com/staylorhead/CADET.git
```

CADET depends opon successful prior installation of:
- [R] (https://www.r-project.org/) and the following packages
    - [lassosum](https://github.com/tshmak/lassosum)
    - [fdrtool](https://cran.r-project.org/web/packages/fdrtool/index.html)
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
