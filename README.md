# GvHD_meta
HSCT GvHD GWAS and gene expression meta-analysis

This repository contains the scripts for running gene-level p-value summarization from Plink-formatted GWAS results using MAGMA, downloading and processing gene expression data from GEO, and rank-based meta-analysis combining these data.

Submitted manuscript: 
Hyvärinen K, Koskela S, Niittyvuopio R, Nihtinen A, Volin L, Salmenniemi U, Putkonen M, Bũno I, Gallardo D, Itälä-Remes M, Partanen J and Ritari J. Meta-analysis of genome-wide association and gene expression studies implicates donor T cell and cytokine pathways in acute GvHD. 

## code (./src)

`main_csc_magma.sh` Shell script for calling subscripts for running MAGMA analysis

`AGVHD.R` R script for running meta-analysis and plotting for aGvHD data
