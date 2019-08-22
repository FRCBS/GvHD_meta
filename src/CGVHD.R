

## Load functions & libs
source('./src/functions.R')


## Data

# pidala   <- fetchGEO('GSE64300')
munagala <- fetchGEO('GSE56495')
hakim    <- fetchGEO('GSE60674')
kohrt    <- fetchGEO('GSE23924')
baron    <- fetchGEO('GSE4624')


## DEGs

munagala.deg <- DEG(munagala, "source_name_ch1", "baseline|disease control", "healthy", T, 
                    gene.id=munagala[['Genes']][, 'Gene symbol'])

hakim.deg    <- DEG(hakim, "disease status:ch1", "CGVHD", "healthy", 
                    gene.id=str_split_fixed(featureData(hakim[[1]])@data$gene_assignment, " // ", 3)[, 2] %>% 
                      str_split_fixed(., " /// ", 3) %>% .[, 1])

kohrt.deg    <- DEG(kohrt, "disease state:ch1", "^cGvHD", "non", gene.id=kohrt[['Genes']][, 'Gene symbol'])

baron.deg    <- DEG(baron, "title", "_cGVHD\\+", "_cGVHD\\-", gene.id=baron[['Genes']][, 'Gene symbol'])


## Rank aggregation with MAGMA GWAS

d.cgvhd.magma.1 <- fread('~/Projects/HSCT_joint/data/magma/IC1_donor_cGvHD12.magma.logistic_ADD_MAGMA.genes.out', data.table=F)
d.cgvhd.magma.2 <- fread('~/Projects/HSCT_joint/data/magma/IC2_3S_donor_cGvHD12.magma.logistic_ADD_MAGMA.genes.out', data.table=F)
d.cgvhd.magma.3 <- fread('~/Projects/HSCT_joint/data/magma/IC3F_donor_cGvHD12.magma.logistic_ADD_MAGMA.genes.out', data.table=F)

d.cgvhd.magma.1 <- d.cgvhd.magma.1 %>% arrange(., P) %>% .$GENE %>% 
  bitr(., fromType='ENTREZID', toType='SYMBOL', OrgDb='org.Hs.eg.db') %>% .$SYMBOL

d.cgvhd.magma.2 <- d.cgvhd.magma.2 %>% arrange(., P) %>% .$GENE %>% 
  bitr(., fromType='ENTREZID', toType='SYMBOL', OrgDb='org.Hs.eg.db') %>% .$SYMBOL

d.cgvhd.magma.3 <- d.cgvhd.magma.3 %>% arrange(., P) %>% .$GENE %>% 
  bitr(., fromType='ENTREZID', toType='SYMBOL', OrgDb='org.Hs.eg.db') %>% .$SYMBOL


# meta donors
d.cgvhd.genelist <- list(d.cgvhd.magma.1, d.cgvhd.magma.2, d.cgvhd.magma.3,
                         #pidala.deg$p.value[pidala.deg$p.value %>% order, ] %>% names %>% emptyFilter, 
                         munagala.deg$p.value[munagala.deg$p.value %>% order, ] %>% names %>% emptyFilter,
                         hakim.deg$p.value[hakim.deg$p.value %>% order, ] %>% names %>% emptyFilter, 
                         kohrt.deg$p.value[kohrt.deg$p.value %>% order, ] %>% names %>% emptyFilter,
                         baron.deg$p.value[baron.deg$p.value %>% order, ] %>% names %>% emptyFilter)
d.cgvhd.genelist.common.genes <- lapply(d.cgvhd.genelist, function(x) x[!(grepl('/', x))])
d.cgvhd.genelist.common <- Reduce(intersect, d.cgvhd.genelist)
d.cgvhd.genelist.common.genes <- lapply(d.cgvhd.genelist, function(x) x[(x %in% d.cgvhd.genelist.common)])

d.cgvhd.genelist.meta <- aggregateRanks(d.cgvhd.genelist, exact=T)
d.cgvhd.genelist.meta$Score %>% hist
d.cgvhd.genelist.meta <- mutate(d.cgvhd.genelist.meta, FDR=p.adjust(Score, method='fdr')) 
d.cgvhd.genelist.meta$Score %>% p.adjust(., method='bonferroni') %>% head
d.cgvhd.genelist.meta %>% head(., 25)

# meta donors gwas and expr separately
d.cgvhd.genelist.meta.gwas <- aggregateRanks(list(d.cgvhd.magma.1, d.cgvhd.magma.2, d.cgvhd.magma.3), exact=T)
d.cgvhd.genelist.meta.expr <- aggregateRanks(list(#pidala.deg$p.value[pidala.deg$p.value %>% order, ] %>% names %>% emptyFilter, 
                                                  munagala.deg$p.value[munagala.deg$p.value %>% order, ] %>% names %>% emptyFilter,
                                                  hakim.deg$p.value[hakim.deg$p.value %>% order, ] %>% names %>% emptyFilter, 
                                                  kohrt.deg$p.value[kohrt.deg$p.value %>% order, ] %>% names %>% emptyFilter,
                                                  baron.deg$p.value[baron.deg$p.value %>% order, ] %>% names %>% emptyFilter), exact=T)
d.cgvhd.genelist.meta.gwas <-  mutate(d.cgvhd.genelist.meta.gwas, FDR=p.adjust(Score, method='fdr')) 
d.cgvhd.genelist.meta.expr <-  mutate(d.cgvhd.genelist.meta.expr, FDR=p.adjust(Score, method='fdr')) 

d.cgvhd.genelist.meta.gwas.bp <- enrichGO(d.cgvhd.genelist.meta.gwas$Name[d.cgvhd.genelist.meta.gwas$FDR<0.1], 
                                     'org.Hs.eg.db', ont='BP', universe=d.cgvhd.genelist.meta.gwas$Name, keyType='SYMBOL')
d.cgvhd.genelist.meta.expr.bp <- enrichGO(d.cgvhd.genelist.meta.expr$Name[d.cgvhd.genelist.meta.expr$FDR<0.1], 
                                          'org.Hs.eg.db', ont='BP', universe=d.cgvhd.genelist.meta.expr$Name, keyType='SYMBOL')


# permutate gene lists

cgvhd.perm <- lapply(1:100, function(iter) {
  print(iter)
  glist.perm <- lapply(d.cgvhd.genelist, function(x) x[sample(1:length(x), length(x), F)])
  perm.meta <- aggregateRanks(glist.perm, exact=T)
  mutate(perm.meta, FDR=p.adjust(Score, method='fdr')) 
})
tmp <- sapply(cgvhd.perm, function(x) x[x$FDR<0.1, 'Name'] %>% as.vector )
d.cgvhd.genelist.meta$Name[d.cgvhd.genelist.meta$FDR<0.1] %>% as.vector %>% length
sapply(tmp, function(x) length(x)>=16) %>% sum



# GO:BP enrichment

d.cgvhd.genelist.meta.bp <- enrichGO(d.cgvhd.genelist.meta$Name[d.cgvhd.genelist.meta$FDR<0.1], 
                                     'org.Hs.eg.db', ont='BP', universe=d.cgvhd.genelist.meta$Name, keyType='SYMBOL')
d.cgvhd.genelist.meta.bp@result$GeneRatio <- d.cgvhd.genelist.meta.bp@result$GeneRatio %>% map(., processGOratios) %>% unlist
d.cgvhd.genelist.meta.bp@result$BgRatio <- d.cgvhd.genelist.meta.bp@result$BgRatio %>% map(., processGOratios) %>% unlist
d.cgvhd.genelist.meta.bp@result <- mutate(d.cgvhd.genelist.meta.bp@result, GeneBgRatio=GeneRatio/BgRatio)
d.cgvhd.genelist.meta.bp@result <- d.cgvhd.genelist.meta.bp@result %>% arrange(., pvalue)
d.cgvhd.genelist.meta.bp@result$ID <- factor(d.cgvhd.genelist.meta.bp@result$ID, levels=d.cgvhd.genelist.meta.bp@result$ID)
d.cgvhd.genelist.meta.bp.flt <- filter(d.cgvhd.genelist.meta.bp@result, p.adjust<0.05)
d.cgvhd.genelist.meta.bp.flt <- d.cgvhd.genelist.meta.bp.flt[!duplicated(d.cgvhd.genelist.meta.bp.flt$geneID), ]

ggplot(d.cgvhd.genelist.meta.bp.flt, aes(pvalue %>% -log10(.), ID, size=GeneBgRatio %>% log)) +
  geom_point()


d.cgvhd.genelist.meta.bp@result %>% filter(., p.adjust<0.05) %>% .[, c('ID', 'pvalue')]

d.cgvhd.genelist.meta.bp@result %>% filter(., p.adjust<0.05) %>% .$Description

write.table(d.cgvhd.genelist.meta.bp@result, './results/D_cGvHD_meta_GO-BP.tsv', sep='\t', row.names=F)
write.table(d.cgvhd.genelist.meta %>% filter(., FDR<0.1), './results/D_cGvHD_meta_FDR015.tsv', sep='\t', row.names=F)


plot((d.cgvhd.genelist.meta.bp@result$GeneRatio %>% map(., processGOratios) %>% unlist /
       d.cgvhd.genelist.meta.bp@result$BgRatio %>% map(., processGOratios) %>% unlist) %>% log,
     d.cgvhd.genelist.meta.bp@result$pvalue %>% -log10(.) )

