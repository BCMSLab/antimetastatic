# load required libraries
library(tidyverse)
library(igraph)

# load data
# knock_tf
knocktf <- read_tsv('data/knock_tf_subset.tsv')
tf_list <- unique(knocktf$TF)

# msg_deg
knockmsg <- read_rds('data/diff_expr_groups.rds')
msg_list <- read_lines('data/msg_list.txt')

# combine
diff_expr <- bind_rows(
    list(
        TF = select(knocktf, KD = TF, ID = Gene, FC = Log2FC, PVAL = P_value),
        MSG = select(knockmsg, KD = group, ID, FC = logFC, PVAL = P.Value, FDR = adj.P.Val)
    ),
    .id = 'Type') %>%
    filter(!KD %in% c('POLR3A', 'KISS1')) %>%
    filter(!ID %in% c('POLR3A', 'KISS1'))

# curated string interactions
string_interactions <- read_csv('data/combined_interactions.csv')

full_edges <- string_interactions %>%
    filter(consensus %in% c(1, -1),
           !grepl('Correlation', Interaction)) %>%
    select(Source.Node = node1,
           Direction = consensus,
           Target.Node = node2)
nrow(full_edges)
# make models list
meta_model <- list()
meta_model$model <- list()

meta_model$model$edges <- as.data.frame(full_edges)

# remove nodes that change with more than one kd condition
meta_model$startNodeDown <- diff_expr %>%
    group_by(KD, ID) %>%
    top_n(n=1, abs(FC)) %>%
    ungroup() %>%
    filter(PVAL < .01 & abs(FC) > .5) %>%
    filter(!ID %in% unique(c(full_edges$Source.Node, full_edges$Target.Node))) %>%
    select(Source.Node = KD,
           Direction = FC,
           Target.Node = ID) %>%
    mutate(Direction = ifelse(Direction > 0, 1, -1),
    Direction * -1) %>%
    as.data.frame() %>%
    with(split(data.frame(nodeLabel = paste0('EXP(', Target.Node, ')'),
                          Direction = Direction),
               Source.Node))

meta_model$g <- meta_model$model$edges %>%
    select(from = Source.Node,
           to = Target.Node,
           weight = Direction) %>%
    graph.data.frame(directed = FALSE)

Hs__CFA__Metastasis__0__0__1 <- meta_model

unlink('NPAModels/data/Hs__CFA__Metastasis__0__0__1.rda')

save(Hs__CFA__Metastasis__0__0__1,
     file = 'NPAModels/data/Hs__CFA__Metastasis__0__0__1.rda')
