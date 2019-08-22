
library(data.table)
library(tidyverse)

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


# list plink assoc files
ls.add <- list.files(pattern='ADD')

# write MAGMA output files
sapply(ls.add, function(FILE) {
	tmp <- fread(FILE, data.table=F)[, c(2, 9)] # SNP P-value
	tmp[, 1] <- str_split_fixed(tmp[, 1], ':', 5)[, 1]
	write.table(tmp, paste0(FILE, '_MAGMA'), quote=F, sep='\t', row.names=F, col.names=c('SNP', 'P'))
})


# find common set of snps
all.snps <- lapply(ls.add, function(FILE) {
	tmp <- fread(FILE, data.table=F)[, 2] # SNP
	str_split_fixed(tmp, ':', 5)[, 1]
})

common.snps <- Reduce(intersect, all.snps)
write.table(common.snps, 'common_snps_rsid.txt', quote=F, col.names=F, row.names=F)


# write common snp assoc files
sapply(ls.add, function(FILE) {
	tmp <- fread(FILE, data.table=F)
	tmp[, 2] <- str_split_fixed(tmp[, 2], ':', 5)[, 1]
	tmp <- tmp[ match(common.snps, tmp[, 2]), ]
	write.table(tmp, paste0(FILE, '_COMMONVAR'), quote=F, sep='\t', row.names=F, col.names=F)
})

# write cleaned snp assoc files
sapply(ls.add, function(FILE) {
	tmp <- fread(FILE, data.table=F)
	tmp[, 2] <- str_split_fixed(tmp[, 2], ':', 5)[, 1]
	tmp <- na.omit(tmp)
	write.table(tmp, paste0(FILE, '_CLEANED'), quote=F, sep='\t', row.names=F, 
		col.names=c('CHR', 'SNP', 'BP', 'A1', 'TEST', 'NMISS', 'BETA', 'STAT', 'P')
})


