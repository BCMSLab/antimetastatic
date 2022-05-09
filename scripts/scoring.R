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

# check no names are duplicated
all(map(diff_expr_groups, ~all(sum(duplicated(toupper(.x))))) == 0)

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

# other cells
fls <- list.files('data/lincs_cells_deg/', full.names = TRUE)
names(fls) <- str_split(fls, '/|\\.', simplify = TRUE)[, 4]

dir.create('data/npa_metastasis/')
dir.create('data/leading_nodes/')

imap(fls, 
     function(x, .y) {
         diff_expr <- read_tsv(x)
         
         diff_expr_groups <- diff_expr %>%
             select(nodeLabel = ID,
                    t,
                    foldChange = logFC) %>%
             as.data.frame() %>%
             with(split(., diff_expr$treatment)) %>%
             map(function(x) {
                 filter(x, !duplicated(x$nodeLabel))
             })
         
         # check no names are duplicated
         all(map(diff_expr_groups, ~all(sum(duplicated(toupper(.x))))) == 0)
         
         # score the metastasis network
         npa_metastasis <- compute_npa(diff_expr_groups,
                                       net.metastasis,
                                       verbose = TRUE)
         
         write_rds(npa_metastasis, paste0('data/npa_metastasis/',.y,'.rds'))
         
         leading_nodes <- as.matrix(npa_metastasis, type = 'leadingnodes') %>%
             melt() %>%
             as_tibble() %>%
             setNames(c('node', 'perturbation', 'value')) %>%
             mutate(signif = grepl('\\*', value)) %>%
             separate(value, into = c('rank', 'direction', 'percentatge'),
                      sep = ' \\! \\(|\\* \\(|  \\(|\\) |\\%',
                      convert = TRUE)
         
         write_rds(leading_nodes, paste0('data/leading_nodes/',.y,'.rds'))
     })