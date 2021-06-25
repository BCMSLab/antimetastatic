# load libraries
library(tidyverse)
library(org.Hs.eg.db)
library(slinky)

# load data
# source: https://docs.google.com/document/d/1q2gciWRhVCAAnlvF2iRLuJ7whrGP6QjpsCMq1yWz7dU/edit
# source: http://bioconductor.org/packages/release/bioc/vignettes/slinky/inst/doc/LINCS-analysis.html#from-the-info-file
sl <- new('Slinky')

# download info
#download(sl, type = 'info')

#files <- GEOquery::getGEOSuppFiles('GSE92742', fetch_files = FALSE)
#write_tsv(files, 'data/GSE92742_files.tsv')

# download.file('https://ftp.ncbi.nlm.nih.gov/geo/series/GSE92nnn/GSE92742/suppl//GSE92742_Broad_LINCS_inst_info.txt.gz',
#               destfile = 'data/GSE92742_Broad_LINCS_inst_info.txt.gz')
# aria2c -x 8 -s 8 https://ftp.ncbi.nlm.nih.gov/geo/series/GSE92nnn/GSE92742/suppl//GSE92742_Broad_LINCS_Level3_INF_mlr12k_n1319138x12328.gctx.gz

user_key <- httr::content(httr::GET("https://api.clue.io/temp_api_key"), 
                          as = "parsed")$user_key

sl <- Slinky(user_key,
             '/home/rstudio/workingon/LINPS/data/GSE92742_Broad_LINCS_Level3_INF_mlr12k_n1319138x12328.gctx',
             '/home/rstudio/workingon/LINPS/data/GSE92742_Broad_LINCS_inst_info.txt.gz')

md <- metadata(sl)

hms <- read_csv('data/hms_gr.csv', guess_max = 1000000)
trt <- hms %>% pull(`Small Molecule Name`) %>% unique()

col.ix <- which(md$cell_id == 'MCF7' & md$pert_iname %in% tolower(trt))
col.ix2 <- which(md$cell_id == 'MCF7' & md$pert_type == 'ctl_vehicle' & md$pert_iname == 'DMSO')

# to be removed later
set.seed(123)
col.ix2 <- sample(col.ix2, 100)

# trt_cp <- readGCTX(sl[, col.ix])
trt_cp <- as(sl[, c(col.ix, col.ix2)], "SummarizedExperiment")

id_symbol <- select(org.Hs.eg.db,
                    rownames(trt_cp),
                    'SYMBOL',
                    'ENTREZID')
all(rownames(trt_cp) == id_symbol$ENTREZID)
rownames(trt_cp) <- id_symbol$SYMBOL

write_rds(trt_cp, 'data/lincs_trt_cp.rds')
