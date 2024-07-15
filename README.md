# mtDNA-CN-ExWAS

Rare variant association testing to mitochondrial copy number in UK Biobank
@Author Vamsee Pillalamarri


This repo contains code (pipeline, scripts, and analysis notebooks) relating to Pillalamarri et al. 2023 HGG Advances "Whole-exome sequencing in 415,422 individuals identifies rare variants associated
with mitochondrial DNA copy number" (PMID 36311265, doi: 10.1016/j.xhgg.2022.100147).


This repository is a snapshot of a private repository with more comprehensive access to a variety of resources for the project along with results from the single- and rare-variant burden testing.


# R Package Requirements

Install BIOCONDUCTOR Install Manager. 

```
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
```

## Aggregate Association Testing  
- GENESIS[https://bioconductor.org/packages/release/bioc/html/GENESIS.html]
- SeqArray[https://bioconductor.org/packages/release/bioc/html/SeqArray.html]
 

Install Bioconductor packages
```
BiocManager::install("SeqArray")
BiocManager::install("GENESIS")
```


## Scripts


## Analysis Notebooks


## Data locations
pVCFs: `/dcl01/arking/data/static/UKBiobank/WES/pVCF`

plink: `/dcl01/arking/data/static/UKBiobank/WES/plink`

note: pVCFs combine multiallelic sites into one record, while plink splits them; the plink *.bim* file reflects this difference

## Gene lists can be found in the "resources" folder in this repo, or on JHPCE at /dcl01/arking/data/static/UKBiobank/gene_lists.

gwas_genes.txt: Genes identified from GWAS
depletion_genes.txt: 15 genes known to be causal for mtDNA depletion syndromes

## Mitocarta gene pathways are located on the server, here: /dcs01/arking/arkinglab/resources/mitoCarta/Human.MitoPathways3.0.gmx

