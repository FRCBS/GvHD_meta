

## functions & libs

library(data.table)
library(tidyverse)
library(GEOquery)
library(qvalue)
library(preprocessCore)
library(RobustRankAggreg)
library(clusterProfiler)
library(biomaRt)
library(ggpubr)
library(locuscomparer); library(cowplot)

source('./src/functions.R')
source('./src/gg_manhattan.R')


## Data

baron     <- fetchGEO('GSE4624')
glauzy    <- fetchGEO('GSE75344')
lupsa     <- fetchGEO('GSE103569')
furlan    <- fetchGEO('GSE73809')
takahashi <- fetchGEO('GSE10572')
buzzeo    <- fetchGEO('GSE7510')



## DEGs

baron.deg2 <- DEG(baron, "title", "_aGVHD\\+", "_aGVHD\\-", gene.id=baron[['Genes']][, 'Gene symbol'])

glauzy.deg <- DEG(glauzy, "characteristics_ch1.1", "grade 2| grade 3", "grade 0", gene.id=glauzy[['Genes']][, 'Gene symbol'])

lupsa.deg <- DEG(lupsa, "characteristics_ch1", "Cutaneous|Gastro", "Absent",
                 gene.id=lupsa[[1]] %>% featureData %>% .@data %>% .$GENE_SYMBOL)

furlan.deg <- DEG(furlan, "description.1", "^GVHD", "^No", 
                  gene.id=str_split_fixed(featureData(furlan[[1]])@data$gene_assignment, " // ", 3)[, 2] %>% 
                    str_split_fixed(., " /// ", 3) %>% .[, 1])

takahashi.deg <- DEG(takahashi, "title", "\\+", "-", gene.id=takahashi[[1]] %>% featureData %>% .@data %>% .$Gene_Symbol)

buzzeo.deg <- DEG(buzzeo, "title", "GVHD", "Ctrl", gene.id=buzzeo[['Genes']][, 'Gene symbol'])


## Rank aggregation together with MAGMA GWAS results

d.agvhd.magma.1 <- fread('~/Projects/HSCT_joint/data/magma/IC1_donor_aGvHD24.magma.logistic_ADD_MAGMA.genes.out', data.table=F)
d.agvhd.magma.2 <- fread('~/Projects/HSCT_joint/data/magma/IC2_3S_donor_aGvHD24.magma.logistic_ADD_MAGMA.genes.out', data.table=F)
d.agvhd.magma.3 <- fread('~/Projects/HSCT_joint/data/magma/IC3F_donor_aGvHD24.magma.logistic_ADD_MAGMA.genes.out', data.table=F)

d.agvhd.magma.1 <- d.agvhd.magma.1 %>% arrange(., P) %>% .$GENE %>% 
  bitr(., fromType='ENTREZID', toType='SYMBOL', OrgDb='org.Hs.eg.db') %>% .$SYMBOL

d.agvhd.magma.2 <- d.agvhd.magma.2 %>% arrange(., P) %>% .$GENE %>% 
  bitr(., fromType='ENTREZID', toType='SYMBOL', OrgDb='org.Hs.eg.db') %>% .$SYMBOL

d.agvhd.magma.3 <- d.agvhd.magma.3 %>% arrange(., P) %>% .$GENE %>% 
  bitr(., fromType='ENTREZID', toType='SYMBOL', OrgDb='org.Hs.eg.db') %>% .$SYMBOL

d.agvhd.genelist <- list(d.agvhd.magma.1, d.agvhd.magma.2, d.agvhd.magma.3,
                         baron.deg2$p.value[baron.deg2$p.value %>% order, ] %>% names %>% emptyFilter, 
                         glauzy.deg$p.value[glauzy.deg$p.value %>% order, ] %>% names %>% emptyFilter,
                         lupsa.deg$p.value[lupsa.deg$p.value %>% order, ] %>% names %>% emptyFilter, 
                         furlan.deg$p.value[furlan.deg$p.value %>% order, ] %>% names %>% emptyFilter,
                         takahashi.deg$p.value[takahashi.deg$p.value %>% order, ] %>% names %>% emptyFilter,
                         buzzeo.deg$p.value[buzzeo.deg$p.value %>% order, ] %>% names %>% emptyFilter)
d.agvhd.genelist <- lapply(d.agvhd.genelist, function(x) x[!(grepl('/', x))])

d.agvhd.genelist.meta <- aggregateRanks(d.agvhd.genelist, exact=T)
d.agvhd.genelist.meta$Score %>% hist
d.agvhd.genelist.meta <- mutate(d.agvhd.genelist.meta, FDR=p.adjust(Score, method='fdr')) 
d.agvhd.genelist.meta$Score %>% p.adjust(., method='bonferroni') %>% head
d.agvhd.genelist.meta %>% head(., 25)
write.table(d.agvhd.genelist.meta %>% filter(FDR<0.1), './results/D_aGvHD_meta_genes.tsv', sep='\t', row.names=F)
  

## GO enrichment
# BG ratio: ontology group genes / all genes

d.agvhd.genelist.meta.bp <- enrichGO(d.agvhd.genelist.meta$Name[d.agvhd.genelist.meta$FDR<0.1], 
                                     'org.Hs.eg.db', ont='BP', universe=d.agvhd.genelist.meta$Name, keyType='SYMBOL')
d.agvhd.genelist.meta.bp@result$GeneRatio <- d.agvhd.genelist.meta.bp@result$GeneRatio %>% map(., processGOratios) %>% unlist
d.agvhd.genelist.meta.bp@result$BgRatio <- d.agvhd.genelist.meta.bp@result$BgRatio %>% map(., processGOratios) %>% unlist
d.agvhd.genelist.meta.bp@result <- mutate(d.agvhd.genelist.meta.bp@result, GeneBgRatio=GeneRatio/BgRatio)
d.agvhd.genelist.meta.bp@result <- d.agvhd.genelist.meta.bp@result %>% arrange(., pvalue)
d.agvhd.genelist.meta.bp@result$ID <- factor(d.agvhd.genelist.meta.bp@result$ID, levels=d.agvhd.genelist.meta.bp@result$ID)

write.table(d.agvhd.genelist.meta.bp@result %>% filter(p.adjust<0.05), './results/D_aGvHD_meta_GO-BP.tsv', sep='\t', row.names=F)

# write full GO result for revigo input
write.table(filter(d.agvhd.genelist.meta.bp@result, p.adjust<0.05)[, c('ID', 'GeneBgRatio')],
            './data/aGvHD_toReviGo_full.tsv', sep='\t', quote=F, row.names=F)


# read manually curated list of enriched immuno-linked BPs
d.agvhd.genelist.meta.bp.flt.imm <- fread('./data/Immunology-focused_GP_BPs.txt', data.table=F)
d.agvhd.genelist.meta.bp.flt.imm <- d.agvhd.genelist.meta.bp.flt.imm %>% arrange(., 1/GeneBgRatio)
d.agvhd.genelist.meta.bp.flt.imm <- mutate(d.agvhd.genelist.meta.bp.flt.imm, Description2=Description %>% str_sub(., 1, 60))
d.agvhd.genelist.meta.bp.flt.imm <- unite(d.agvhd.genelist.meta.bp.flt.imm, 'BP', Description2, ID, sep=' ', remove=F)
d.agvhd.genelist.meta.bp.flt.imm$BP <- factor(d.agvhd.genelist.meta.bp.flt.imm$BP, levels=d.agvhd.genelist.meta.bp.flt.imm$BP %>% rev)

p3.2 <- ggplot(d.agvhd.genelist.meta.bp.flt.imm[1:30, ], aes(GeneBgRatio, BP, size=p.adjust %>% log10 %>% `*`(-1))) +
  geom_point(alpha=0.5, stroke=0.2, color='red') + 
  xlab('Fold enrichment') + ylab('GO:BP') +
  guides(size=guide_legend(title=expression('-log'[10]*'(FDR)'))) +
  theme_minimal() +
  coord_cartesian(clip="off") +
  scale_size(range=c(2, 6)) +
  ggtitle('') +
  theme(legend.position='bottom', axis.text.y=element_text(size=8.1), 
        panel.grid.minor.x=element_blank(),
        panel.grid.minor.y=element_blank(),
        panel.grid.major.x=element_blank(),
        panel.grid.major = element_line(size=0.2),
        axis.line=element_line(colour="black", size=.2))


## Meta donors gwas and expr separately

d.agvhd.genelist.meta.gwas <- aggregateRanks(list(d.agvhd.magma.1, d.agvhd.magma.2, d.agvhd.magma.3), exact=T)
d.agvhd.genelist.meta.expr <- aggregateRanks(list(baron.deg2$p.value[baron.deg2$p.value %>% order, ] %>% names %>% emptyFilter, 
                                                  glauzy.deg$p.value[glauzy.deg$p.value %>% order, ] %>% names %>% emptyFilter,
                                                  lupsa.deg$p.value[lupsa.deg$p.value %>% order, ] %>% names %>% emptyFilter, 
                                                  furlan.deg$p.value[furlan.deg$p.value %>% order, ] %>% names %>% emptyFilter,
                                                  takahashi.deg$p.value[takahashi.deg$p.value %>% order, ] %>% names %>% emptyFilter,
                                                  buzzeo.deg$p.value[buzzeo.deg$p.value %>% order, ] %>% names %>% emptyFilter), exact=T)
d.agvhd.genelist.meta.gwas <-  mutate(d.agvhd.genelist.meta.gwas, FDR=p.adjust(Score, method='fdr')) 
d.agvhd.genelist.meta.expr <-  mutate(d.agvhd.genelist.meta.expr, FDR=p.adjust(Score, method='fdr')) 

d.agvhd.genelist.meta.gwas.bp <- enrichGO(d.agvhd.genelist.meta.gwas$Name[d.agvhd.genelist.meta.gwas$FDR<0.1], 
                                          'org.Hs.eg.db', ont='BP', universe=d.agvhd.genelist.meta.gwas$Name, keyType='SYMBOL')
d.agvhd.genelist.meta.expr.bp <- enrichGO(d.agvhd.genelist.meta.expr$Name[d.agvhd.genelist.meta.expr$FDR<0.1], 
                                          'org.Hs.eg.db', ont='BP', universe=d.agvhd.genelist.meta.expr$Name, keyType='SYMBOL')

d.agvhd.comp <- data.frame(SignifGenes=c(d.agvhd.genelist.meta.gwas %>% filter(., FDR<0.1) %>% nrow, 
                                         d.agvhd.genelist.meta.expr %>% filter(., FDR<0.1) %>% nrow,
                                         d.agvhd.genelist.meta %>% filter(., FDR<0.1) %>% nrow),
                           UniqueBPs=c(0, 
                                       d.agvhd.genelist.meta.expr.bp@result %>% filter(., p.adjust<0.05) %>% 
                                         .$geneID %>% unique %>% length,
                                       d.agvhd.genelist.meta.bp@result %>% filter(., p.adjust<0.05) %>% 
                                         .$geneID %>% unique %>% length),
                           Data=c('GWAS', 'Expression', 'GWAS + \nExpression'))

p1.1 <- ggplot(d.agvhd.comp %>% melt %>% filter(., variable=='SignifGenes'), aes(Data, value)) +
  geom_bar(stat='identity', fill='grey60', alpha=0.8, width=0.4) +
  xlab('') + ylab('# Genes (FDR<0.1)') +
  ggtitle('') +
  theme_minimal() + theme(plot.title=element_text(size=10)) +
  scale_y_continuous(expand = c(0, 0)) + 
  theme(axis.line = element_line(colour="black", size=.2),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.text.x = element_text(color='black')) 

p1.2 <- ggplot(d.agvhd.comp %>% melt %>% filter(., variable=='UniqueBPs'), aes(Data, value)) +
  xlab('') + ylab('# GO:BPs (FDR<0.05)') +
  geom_bar(stat='identity', fill='grey60', alpha=0.8, width=0.4) +
  ggtitle('') +
  theme_minimal() + theme(plot.title=element_text(size=10)) +
  scale_y_continuous(expand = c(0, 0)) + 
  theme(axis.line = element_line(colour="black", size=.2),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.text.x = element_text(color='black')) 



## Permutation analysis of gene lists

agvhd.perm <- lapply(1:100, function(iter) {
  print(iter)
  glist.perm <- lapply(d.agvhd.genelist, function(x) x[sample(1:length(x), length(x), F)])
  perm.meta <- aggregateRanks(glist.perm, exact=T)
  mutate(perm.meta, FDR=p.adjust(Score, method='fdr')) 
})
agvhd.perm.flt <- sapply(agvhd.perm, function(x) x[x$FDR<0.1, 'Name'] %>% as.vector)
agvhd.perm.flt.len <- sapply(agvhd.perm.flt, length)
agvhd.perm.overp <- sapply(agvhd.perm, function(x) sum(x$Score<0.0001392108))

# histogram of the number of fdr<0.1 genes
p2 <- qplot(agvhd.perm.flt.len, alpha=0.8, binwidth=2) + 
  geom_line(aes(x=c(51, 51), y=c(0, 30)), color='red', linetype='solid') + 
  xlab('# Genes (FDR<0.1)') + ylab('# permutations') +
  theme_minimal() +
  ggtitle('') +
  theme(legend.position="none", plot.title=element_text(size=10)) +
  scale_y_continuous(expand=c(0, 0)) + 
  scale_colour_manual(values='grey60') +
  theme(axis.line = element_line(colour="black", size=.2),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) 


# GO:BP enrichment on permutated genelists
agvhd.perm.bp <- lapply(agvhd.perm.flt, function(i) {
  print(i)
  if(length(i) > 0) {
    enrichGO(i, 'org.Hs.eg.db', ont='BP', universe=d.agvhd.genelist.meta$Name, keyType='SYMBOL')
  } else 0
})
agvhd.perm.bp.num <- lapply(agvhd.perm.bp, function(i) {
  tryCatch(i@result %>% filter(p.adjust<0.05) %>% .$geneID %>% unique %>% length, error=function(e) 0) 
})

d.agvhd.genelist.meta.bp@result %>% filter(p.adjust<0.05) %>% .$pvalue %>% max
agvhd.perm.bp.overp <- sapply(agvhd.perm.bp, function(x) 
  tryCatch( (x@result[x@result$pvalue<0.006885418, 'geneID'] %>% unique %>% length), error=function(e) 0))

agvhd.perm.bp.num.mean <- sapply(1:1000, function(iter) {
  tt <- agvhd.perm.bp.num %>% unlist
  tt[sample(1:length(tt), length(tt), replace=T)] %>% mean
}) 
agvhd.perm.bp.num.mean %>% mean


# histogram of the number of fdr<0.05 BPs
p4 <- qplot(agvhd.perm.bp.num %>% unlist, alpha=0.8, binwidth=4) + 
  geom_line(aes(x=c(156, 156), y=c(0, 50)), color='red', linetype='solid') + 
  xlab('# GO:BPs (FDR<0.05)') + ylab('# permutations') +
  ggtitle('') +
  theme_minimal() +
  theme(legend.position="none", plot.title=element_text(size=10)) +
  scale_y_continuous(expand=c(0, 0)) + 
  scale_colour_manual(values='grey60') +
  theme(axis.line = element_line(colour="black", size=.2),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) 

# plot
cairo_pdf('./results/aGvHD_meta_comparison.pdf', height=5, width=5)
ggarrange(
  ggarrange(p1.1, p1.2, ncol=2),
  ggarrange(p2, p4, ncol=2),
  nrow=2, ncol=1, labels=c('a', 'b'), align='v')
dev.off()


## Manhattan plot
source('./src/gg_manhattan.R')
# gene coordinates from biomart GRCh37
ensembl <- useMart(biomart="ENSEMBL_MART_ENSEMBL", host="grch37.ensembl.org", path="/biomart/martservice",
                   dataset="hsapiens_gene_ensembl")
d.agvhd.genelist.meta.ensembl <- getBM(attributes=c('wikigene_name', 'wikigene_description', 'ensembl_gene_id', 
                                                    'chromosome_name','start_position', 'end_position') , 
      filters='wikigene_name', values=d.agvhd.genelist.meta$Name, mart=ensembl)

d.agvhd.genelist.meta.ensembl <- d.agvhd.genelist.meta.ensembl[str_length(d.agvhd.genelist.meta.ensembl$chromosome_name) < 5, ]
d.agvhd.genelist.meta.ensembl <- inner_join(d.agvhd.genelist.meta, d.agvhd.genelist.meta.ensembl, by=c('Name'='wikigene_name'))
d.agvhd.genelist.meta.manh <- d.agvhd.genelist.meta.ensembl[, c('Name', 'chromosome_name', 'start_position', 'Score')]
colnames(d.agvhd.genelist.meta.manh) <- c('SNP', 'CHR', 'BP', 'P')
d.agvhd.genelist.meta.manh <- d.agvhd.genelist.meta.manh[!(d.agvhd.genelist.meta.manh$CHR %>% grepl("Y|MT|X", .)), ]
d.agvhd.genelist.meta.manh$CHR <- d.agvhd.genelist.meta.manh$CHR %>% factor(., levels=c(1:22))

cairo_pdf('./results/aGvHD_meta_manhattan.pdf', height=7, width=18)
gg.manhattan(df=d.agvhd.genelist.meta.manh, hlight=NA, sugg=NA,
             threshold=d.agvhd.genelist.meta[which(d.agvhd.genelist.meta$FDR<0.1)+1, 'Score'] %>% max, 
             sig=d.agvhd.genelist.meta[which(d.agvhd.genelist.meta$FDR<0.1), 'Score'] %>% max, 
             col=c('#386cb0', 'skyblue1'), title='', pointsize=3.0, ylims=c(0, 6.9))
dev.off()


## ReviGO output plot
source('./src/aGvHD_REVIGO.r')
one.data <- inner_join(one.data, d.agvhd.genelist.meta.bp@result[, c('ID', 'pvalue', 'Count')], by=c('term_ID'='ID'))
one.data <- dplyr::rename(one.data, Uniqueness=uniqueness)
write.table(one.data, './data/Revigo_plot_items.tsv', sep='\t', row.names=F)
ex <- one.data[!grepl("brain|neural|neuro|glia|glio|dendrit|axon|muscle|ossi|ovulation|osteo|chondro|nephro", one.data$description), ]
ex <- ex[grepl(
"lympho|leuk|immu|kappa|antibiotic|bacteria|viral|virus|peptidoglycan|cytokine|remodel|differentiation|inter|toll|junction|steroid|vitamin", 
               ex$description), ]
ex <- ex[!(ex$term_ID %in% 
             c('GO:0043124', 'GO:0034105', 'GO:1904996', 'GO:0070106', 'GO:0035723', 'GO:0071104', 'GO:0070757', 'GO:0070672')), ]

p5 <- ggplot(one.data, aes(plot_X, plot_Y, size=Count, fill=Uniqueness)) +
  geom_point(shape=21, alpha=0.6, color=I(alpha("black", 0.85)), stroke=0.1) +
  scale_fill_gradientn(colors=c("#386cb0", "green", "yellow", "red"), 
                       breaks=c(0.75, 0.85, 0.95), labels=c(0.75, 0.85, 0.95)) +
  scale_size(range=c(5, 15), guide='none') +
  geom_label_repel(data=ex, aes(plot_X, plot_Y, label=description), color=I(alpha("black", 1)), fill='white', 
                   size=1.95, alpha=0.75, force=4, segment.alpha=0.5, segment.size=0.2) +
  xlab('Semantic similarity dimension 1') + ylab('Semantic similarity dimension 2') +
  ggtitle('') +
  theme_minimal() +
  coord_cartesian(clip="off") +
  theme(axis.line = element_line(colour="black", size=.2),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        legend.position='bottom') 
p5

cairo_pdf('./results/agvhd_GOs.pdf', height=5, width=10)
ggarrange(p3.2, p5, ncol=2, align='h', labels=c('a', 'b'), widths=c(1, 0.8))
dev.off()



## Expression effect sizes and directions

baron.expr <- getExpr(baron, (d.agvhd.genelist.meta %>% filter(FDR<0.1) %>% .$Name), 
                      'title', "_aGVHD\\+", "_aGVHD\\-",
                      gene.id=baron[['Genes']][, 'Gene symbol'])

glauzy.expr <- getExpr(glauzy, (d.agvhd.genelist.meta %>% filter(FDR<0.1) %>% .$Name), 
                       'characteristics_ch1.1', "grade 2| grade 3", "grade 0",
                       gene.id=glauzy[['Genes']][, 'Gene symbol'])

lupsa.expr <- getExpr(lupsa, (d.agvhd.genelist.meta %>% filter(FDR<0.1) %>% .$Name), 
                      "characteristics_ch1", "Cutaneous|Gastro", "Absent",
                      gene.id=lupsa[[1]] %>% featureData %>% .@data %>% .$GENE_SYMBOL)

furlan.expr <- getExpr(furlan, (d.agvhd.genelist.meta %>% filter(FDR<0.1) %>% .$Name), 
                       "description.1", "^GVHD", "^No",
                       gene.id=str_split_fixed(featureData(furlan[[1]])@data$gene_assignment, " // ", 3)[, 2] %>% 
                         str_split_fixed(., " /// ", 3) %>% .[, 1])

takahashi.expr <- getExpr(takahashi, (d.agvhd.genelist.meta %>% filter(FDR<0.1) %>% .$Name), 
                          "title", "\\+", "-",
                          gene.id=takahashi[[1]] %>% featureData %>% .@data %>% .$Gene_Symbol)

buzzeo.expr <- getExpr(buzzeo, (d.agvhd.genelist.meta %>% filter(FDR<0.1) %>% .$Name), 
                       "title", "GVHD", "Ctrl",
                       gene.id=buzzeo[['Genes']][, 'Gene symbol'])


expr.genes.cases <- (d.agvhd.genelist.meta %>% filter(FDR<0.1) %>% .$Name) %>% map(function(gene) {
  list(baron.expr, glauzy.expr, lupsa.expr, furlan.expr, takahashi.expr, buzzeo.expr) %>% map(function(x) {
    scale((x$cases[(x$cases %>% rownames) %in% gene, ] / mean(x$controls[(x$controls %>% rownames) %in% gene, ], na.rm=T)))
  }) %>% unlist %>% na.omit
})
names(expr.genes.cases) <- (d.agvhd.genelist.meta %>% filter(FDR<0.1) %>% .$Name)
expr.genes.cases <- expr.genes.cases %>% melt
expr.genes.cases$L1 <- factor(expr.genes.cases$L1, ordered=T,
                              levels=group_by(expr.genes.cases, L1) %>% summarise(., mean=median(value)) %>% arrange(mean) %>% .$L1)

cairo_pdf('./results/D_aGvHD_expr_meta_case_effects.pdf', width=5, height=6.5)
ggplot(expr.genes.cases, aes(L1, value)) +
  geom_boxplot(fill='grey', outlier.shape=1, outlier.size=0.6, size=0.1, alpha=0.4, width=0.6) +
  ylim(c(-2,2)) +
  geom_hline(yintercept=0, color='red', alpha=0.8, size=0.4) +
  stat_summary(geom="crossbar", width=0.6, fatten=2, color="grey40", 
               fun.data=function(x){ return(c(y=median(x), ymin=median(x), ymax=median(x))) }) +
  xlab('Gene') + ylab('SD') +
  theme_minimal() +
  theme(axis.line = element_line(colour="black", size=.2),
        axis.text.y=element_text(size=8.4),
        panel.grid = element_blank(),
        panel.border = element_blank(), 
        panel.background = element_blank()) +
  coord_flip() 
dev.off()


