# load libraries
library(SummarizedExperiment)
library(limma)
library(tidyverse)

# load data
trt_cp <- read_rds('data/lincs_trt_cp.rds')
trt_cp$pert_iname <- relevel(factor(trt_cp$pert_iname), ref = 'DMSO')

# deg
mod <- model.matrix(~pert_iname, colData(trt_cp))
colnames(mod) <- levels(trt_cp$pert_iname)

fit <- lmFit(assay(trt_cp), mod)
fit <- eBayes(fit)

trt <- colnames(fit$coefficients)[-1]
names(trt) <- trt

res <- map_df(trt, function(x) {
    topTable(fit, x, number = Inf)
}, .id = 'treatment') %>%
    as_tibble()

write_tsv(res, 'data/lincs_trt_cp_deg.tsv')

# other cells
fls <- list.files('data/lincs_cells/', full.names = TRUE)
names(fls) <- str_split(fls, '/|\\.', simplify = TRUE)[, 4]

dir.create('data/lincs_cells_deg')

imap(fls,
    function(x, .y) {
      trt_cp <- read_rds(x)
      trt_cp$pert_iname <- relevel(as.factor(trt_cp$pert_iname), ref = 'DMSO')
      
      mod <- model.matrix(~pert_iname, colData(trt_cp))
      colnames(mod) <- levels(trt_cp$pert_iname)

      fit <- lmFit(assay(trt_cp), mod)
      fit <- eBayes(fit)
      
      trt <- colnames(fit$coefficients)[-1]
      names(trt) <- trt
      
      res <- map_df(trt, function(x) {
        topTable(fit, x, number = Inf)
      }, .id = 'treatment') %>%
        as_tibble()
      
      write_tsv(res, paste0('data/lincs_cells_deg/', .y, '.tsv'))
    }
)
