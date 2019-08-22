
## CSC

## R; generate magma input files
cd /proj/2000066/Temp/HSCT_logreg_JointGWAS/assoc

# plink logistic output
#     CHR       Chromosome
#     SNP       SNP identifier
#     BP        Physical position (base-pair)
#     A1        Tested allele (minor allele by default) 
#     TEST      Code for the test (see below)
#     NMISS     Number of non-missing individuals included in analysis
#     BETA/OR   Regression coefficient (--linear) or odds ratio (--logistic)
#     STAT      Coefficient t-statistic 
#     P         Asymptotic p-value for t-statistic

module load r-env
Rscript process_PLINK_logistic.R

## MAGMA runs
cd /proj/2000066/Temp/HSCT_logreg_JointGWAS/magma
magma_run.sh

## R; common genes in magma analyses
process_MAGMA_output.R

