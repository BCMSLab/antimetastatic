# load libraries
library(tidyverse)
library(Biobase)
library(limma)
library(sva)

# load data
msg <- read_rds('data/msg_kd_eset.rds')
ids <- read_lines('data/msg_list.txt')

tt <- map_df(ids, function(x) {
    e <- (msg)[, msg$group %in% c('None', x)]

    mat <- exprs(e)
    mat <- mat[rowSums(mat) > 10,]
    mat <- log2(mat)

    e$group <- droplevels(e$group)
    mod <- model.matrix(~e$group)

    fit <- lmFit(mat, mod)
    fit <- eBayes(fit)

    topTable(fit, number = Inf, genelist = rownames(mat)) %>%
            as_tibble() %>%
            mutate(group = x)
})

write_rds(tt, 'data/diff_expr_groups.rds')
