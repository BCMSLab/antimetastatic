library(tidyverse)
library(reshape2)
library(NPA)

string_interactions <- read_csv('data/combined_interactions.csv') %>%
    filter(consensus %in% c(-1, 1),
           !grepl('Correlation', Interaction))

ll1 <- string_interactions %>%
    group_by(node1) %>%
    summarise(node2 = list(node2),
              expect = list(consensus)) %>%
    group_split(node1)

# npa
npa_metastasis <- read_rds('data/npa_metastasis.rds')

ll2 <- coefficients(npa_metastasis, type = 'nodes') %>%
    melt() %>%
    setNames(c('node', 'drug', 'coeff')) %>%
    # mutate(coeff = ifelse(coeff > 0, 1, -1)) %>%
    group_split(drug) 

df <- ll2 %>%
    map_df(function(x) {
        map_df(ll1, function(y) {
            upstream_coef <- x$coeff[x$node == y$node1]
            upstream <- ifelse(upstream_coef > 0, 1, -1)
            
            downstream_coef <- x$coeff[x$node %in% unlist(y$node2)]
            downstream <- ifelse(downstream_coef > 0, 1, -1)
            tibble(y) %>%
                mutate(upstream = upstream,
                       upstream_coef = upstream_coef,
                       real = list(downstream),
                       downstream_coef = list(downstream_coef),
                       drug = unique(x$drug))
        })
    }) 

df <- df[lengths(df$expect) == lengths(df$real),]
df$con <- pmap_dbl(list(df$expect, df$real, df$upstream),
                function(x, y, z) mean(x == (y * z)))

write_rds(df, 'data/true_labels.rds')
