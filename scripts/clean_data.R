# load required libraries
library(GEOquery)
library(tidyverse)
library(WGCNA)
allowWGCNAThreads(4)

# load data
eset <- getGEO('GSE53668', destdir = 'data/')[[1]]
md <- read_csv('data/md.csv')
md <- as.data.frame(md)
rownames(md) <- colnames(eset)

mat <- exprs(eset)
fd <- fData(eset) %>% as.data.frame()

mat_collapsed <- collapseRows(mat,
                              fd$`Gene Symbol`,
                              rownames(mat))
ind <- grepl('///', rownames(mat_collapsed$datETcollapsed))

mat2 <- as.matrix(mat_collapsed$datETcollapsed[!ind,])

mds <- cmdscale(dist(t(log2(mat2 + 1))))

e2 <- ExpressionSet(mat2,
                    phenoData = new('AnnotatedDataFrame', as.data.frame(md)))

e2$target <- relevel(as.factor(e2$target), ref = 'Scramble')
# e2$time <- relevel(as.factor(e2$time), ref = '0')
e2$time <- as.integer(as.character(e2$time))
e2$type <- relevel(as.factor(e2$type), ref = 'WT')
e2$group <- relevel(as.factor(e2$group), ref = 'None')

write_rds(e2, 'data/msg_kd_eset.rds')
