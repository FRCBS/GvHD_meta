
library(data.table)
library(tidyverse)

# list magma output files
ls.mag <- list.files(pattern='genes.out')

# find common set of genes
all.genes <- lapply(ls.mag, function(FILE) fread(FILE, data.table=F)$GENE )
common.genes <- Reduce(intersect, all.genes)
write.table(common.genes, 'common_genes_magma.txt', quote=F, col.names=F, row.names=F)

# generate common gene magma output files
sapply(ls.mag, function(FILE) {
	tmp <- fread(FILE, data.table=F)
	tmp <- tmp[ match(common.genes, tmp$GENE), ]
	write.table(tmp, paste0(FILE, '_COMMONGENE'), quote=F, sep=' ', row.names=F)
})


