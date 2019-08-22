#!/bin/bash

for FILE in `ls -1 /proj/2000066/Temp/HSCT_logreg_JointGWAS/assoc/*MAGMA`; do
  OUTFILE=`echo $FILE | sed -e "s/assoc/magma/g"`
  /proj/2000066/tools/MAGMA/bin/magma --bfile /proj/2000066/tools/MAGMA/data/g1000_eur --pval $FILE N=1000 --gene-annot /proj/2000066/tools/MAGMA/data/g1000_eur_37_annotation.genes.annot --out $OUTFILE

done

