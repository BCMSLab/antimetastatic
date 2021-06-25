# load libraries
library(tidyverse)

# download data
url <- 'http://www.licpathway.net/KnockTF/download_FC/differential%20expression%20of%20genes%20in%20all%20datasets.txt'
file <- 'data/differential expression of genes in all datasets.txt'

if (!file.exists(file)) {
    download.file(url, file)
}

tf_list <- read_lines('data/tf_list.txt')

knocktf <- read_tsv(file)

url2 <- 'http://www.licpathway.net/KnockTF/search/search_biosample_result.php?sample_tissue_type=Mammary_gland&sample_biosample_name=MCF7&sample_biosample_type=Cell+line&sample_source=All#'
file2 <- 'data/KnockTF-Search result of Tissue Type.csv'

if (!file.exists(file2)) {
    download.file(url2, file2)
}

md <- read_csv(file2)
dataset_ids <- md %>%
    filter(TF %in% tf_list, `Biosample Name` == 'MCF7') %>%
    pull(`Dataset ID`)

dataset_ids <- c("DataSet_01_51", 
                 "DataSet_01_64","DataSet_01_86",
                 "DataSet_01_99", "DataSet_01_111",
                 "DataSet_01_120", "DataSet_01_139",
                 "DataSet_01_150","DataSet_01_182","DataSet_01_189",
                 "DataSet_01_223","DataSet_02_206")
all(dataset_ids %in% unique(knocktf$Sample_ID))

ind <- which(knocktf$Sample_ID %in% dataset_ids)
knocktf[ind,] %>%
    write_tsv('data/knock_tf_subset.tsv')


