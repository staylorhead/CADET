# CADET: Enhanced transcriptome-wide association analyses in admixed samples using eQTL summary data 

## Introduction

CADET is a tool that enables powerful transcriptome-wide association study (TWAS) of admixed cohorts leveraging the local-ancestry (LA) information of the cohort along with summary-level eQTL data from reference panels of different ancestral groups. CADET combines multiple polygenic risk score models with the summary-level eQTL reference data to predict LA-aware genetically-regulated gene expression in target admixed samples. This gene expression can then be used in TWAS to test for gene-level association with a given phenotype.

For questions or issues related to this software, please contact Taylor Head (<sthead@mdanderson.org>).

## Installation 

CADET depends on [OTTERS](https://www.nature.com/articles/s41467-023-36862-w). Please follow the installation instructions [here](https://github.com/daiqile96/OTTERS) and create the OTTERS environment.

```bash
git clone https://github.com/staylorhead/CADET.git
```

CADET also depends opon installation of 
- [BGZIP](http://www.htslib.org/doc/bgzip.html)
- [TABIX](http://www.htslib.org/doc/tabix.html)
- [PLINK 1.9](https://www.cog-genomics.org/plink/)
- Python modules/libraries
    - [pandas 1.4.4](https://pandas.pydata.org)
    - [scipy 1.7.3](https://scipy.org)
    - [numpy 1.21.5](https://numpy.org)
    - [pysam 0.19.1](https://pysam.readthedocs.io/en/latest/api.html) 

### Input Files

### Example Code
