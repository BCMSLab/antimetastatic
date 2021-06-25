# loading required libraries
library(tidyverse)
library(reshape2)

devtools::install('NPAModels/')

library(NPA)
library(NPAModels)
preprocessNetworks()

net.metastasis <- load_model('Hs', 'CFA', 'Metastasis')

# load data
diff_expr <- read_tsv('data/lincs_trt_cp_deg.tsv')

diff_expr_groups <- diff_expr %>%
    select(nodeLabel = ID,
           t,
           foldChange = logFC) %>%
    as.data.frame() %>%
    with(split(., diff_expr$treatment)) %>%
    map(function(x) {
        filter(x, !duplicated(x$nodeLabel))
    })

# ind <- Reduce(intersect, map(diff_expr_groups, ~.x$nodeLabel))
# 
# diff_expr_groups <- map(diff_expr_groups, function(x) {
#     df <- x
#     df <- df[df$nodeLabel %in% ind,]
#     df <- df[!duplicated(toupper(df$nodeLabel)),]
#     rownames(df) <- df$nodeLabel
#     df
# })

# check no names are duplicated
all(map(diff_expr_groups, ~all(sum(duplicated(toupper(.x))))) == 0)

# map(models, function(x) {
#     nn <- x$get_data()$model[[1]] %>% nrow()
#     nt <- x$get_data()$startNodeDown %>% map(nrow) %>%  unlist() %>% sum()
#     sum(nn, nt)
# })

# # score multiple models
# npa <- compute_npa_list(diff_expr_groups,
#                         models[c(2, 4, 6)],
#                         verbose = TRUE)
# write_rds(npa, 'data/npa.rds')
# 
# # calculate impact factor
# bif <- get_bif(npa)
# write_rds(bif, 'data/bif.rds')

# score the metastasis network
npa_metastasis <- compute_npa(diff_expr_groups,
                              net.metastasis,
                              verbose = TRUE)

write_rds(npa_metastasis, 'data/npa_metastasis.rds')

leading_nodes <- as.matrix(npa_metastasis, type = 'leadingnodes') %>%
    melt() %>%
    as_tibble() %>%
    setNames(c('node', 'perturbation', 'value')) %>%
    mutate(signif = grepl('\\*', value)) %>%
    separate(value, into = c('rank', 'direction', 'percentatge'),
             sep = ' \\! \\(|\\* \\(|  \\(|\\) |\\%',
             convert = TRUE)

write_rds(leading_nodes, 'data/leading_nodes.rds')
