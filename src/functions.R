

## Libs

library(data.table)
library(tidyverse)
library(GEOquery)
library(qvalue)
library(preprocessCore)
library(RobustRankAggreg)
library(clusterProfiler)


## Functions

# Remove empty entries
emptyFilter <- function(x) x[!(x=='')]

# Getting GEO data
fetchGEO <- function(x) { 
  tmp <- getGEO(x, AnnotGPL=T, destdir='./data', GSEMatrix=T)[[1]]
  out <- list(tmp,
              fData(tmp)[, 3:4],
              exprs(tmp),
              tmp %>% phenoData %>% .@data)
  names(out) <- c('Data', 'Genes', 'Expression', 'Samples')
  out
}

# GO enrichment; ratios
processGOratios <- function(x) {
  tmp <- str_split(x, '/')[[1]]
  as.numeric(tmp[1]) / as.numeric(tmp[2])
}

# Differentially expressed genes
DEG <- function(d, target, target.dis, target.ctr, log.transform=F, gene.id=NULL) {
    
  if(log.transform==F) tmp <- d[['Expression']] else tmp <- log2(d[['Expression']]) 
  
  print(length(gene.id)) ; print(dim(tmp))
  tmp <- tapply(1:nrow(tmp), gene.id, function(x) colMeans(tmp[x, ] %>% data.frame, na.rm=T) ) %>% do.call(rbind, .)
  
  pheno <- rep(NA, nrow(d[['Samples']]))
  pheno[grepl(target.dis, d[['Samples']][, target])] <- 1
  pheno[grepl(target.ctr, d[['Samples']][, target])] <- 0
  pheno[d[['Samples']][, target] %>% is.na] <- NA
  
  tmp <- tmp[, !is.na(pheno)]
  pheno <- pheno[!is.na(pheno)]
  print(pheno)  

  design <- model.matrix(~0+as.factor(pheno))
  colnames(design) <- c("GvHD","control")
  contrast.matrix <- makeContrasts(GvHD_v_control=GvHD-control, levels=design)
  
  fit <- lmFit(tmp, design)
  fit2 <- contrasts.fit(fit, contrast.matrix)
  eBayes(fit2)
  
}


# Extract expression values
getExpr <- function(d, genes, sam.header, controls, cases, gene.id) {
  gene.ind <- gene.id %in% genes
  express <- d[['Expression']][gene.ind, ]
  rownames(express) <- gene.id[gene.ind]
  if(any(duplicated(rownames(express)))) {
    express <- tapply(1:nrow(express), rownames(express), function(x)
      colMeans(express[x, ] %>% data.frame, na.rm=T) ) %>% do.call(rbind, .)
  } else express <- express
  control.ind <- grepl(controls, d[['Samples']][, sam.header])
  cases.ind <- grepl(cases, d[['Samples']][, sam.header])
  out <- list(express[, control.ind], express[, cases.ind]) 
  names(out) <- c('controls', 'cases')
  out
}

