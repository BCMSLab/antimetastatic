library(tidyverse)

# source: https://lincs.hms.harvard.edu/db/datasets/20268
hms1 <- read_csv('data/dataset_20268_20210412113444.csv')

# source: https://lincs.hms.harvard.edu/db/datasets/20269
hms2 <- read_csv('data/dataset_20269_20210412114601.csv')

hms <- full_join(hms1, hms2)
write_csv(hms, 'data/hms_gr.csv')


hms <- read_csv('data/hms_gr.csv')
