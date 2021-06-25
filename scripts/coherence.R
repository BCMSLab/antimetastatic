library(tidyverse)
library(reshape2)
library(NPA)
library(igraph)

string_interactions <- read_csv('data/combined_interactions.csv') %>%
    filter(consensus %in% c(-1, 1),
           !grepl('Correlation', Interaction))

g <- graph_from_data_frame(
    tibble(
        from = string_interactions$node1,
        to = string_interactions$node2,
        weight = string_interactions$consensus
    ),
    directed = TRUE
)

paths <- all_simple_paths(g, from = 'PEBP1', mode = 'in')
names(paths) <- 1:length(paths)

paths2 <- list()
for (i in seq_along(paths)) {
    nms <- names(paths[[i]])
    nms <- nms[length(nms):1]
    e <- list()
    j <- 0
    while(j < length(nms)-1) {
        j <- j + 1
        e[[j]] <- (c(nms[j], nms[j+1]))
    }
    df <- as.data.frame(do.call('rbind', e))
    names(df) <- paste0('node', 1:2)
    paths2[[i]] <- df
}
names(paths2) <- 1:length(paths2)
paths2 <- bind_rows(paths2, .id = 'path') %>% as_tibble()

# npa
npa_metastasis <- read_rds('data/npa_metastasis.rds')

coeffs <- coefficients(npa_metastasis, type = 'nodes') %>%
    melt() %>%
    setNames(c('node', 'drug', 'coeff'))

coherence <- left_join(paths2, select(string_interactions, node1, node2, consensus)) %>%
    left_join(coeffs, by = c('node1'='node')) %>%
    left_join(coeffs, by = c('node2'='node', 'drug' = 'drug')) %>%
    mutate(node1_perturb = coeff.x, node2_perturb = coeff.y) %>%
    mutate_at(vars(starts_with('coeff')), function(x) ifelse(x > 0, 1, -1)) %>%
    group_split(path, drug) %>%
    map_df(function(x) {
        v <- c(x$coeff.y[1], x$consensus)
        r <- vector()
        i <- 1
        while (i < length(v)) {
            r[i] <- v[i] * v[i+1]
            i <- i +1
        }
        x$observed <- r
        x$coherence <- as.numeric(x$coeff.y == x$observed)
        x
    })

write_csv(coherence, 'data/coherence.csv')
