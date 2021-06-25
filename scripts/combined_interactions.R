# tf deg
knocktf <- read_tsv('data/knock_tf_subset.tsv')
tf_list <- unique(knocktf$TF)

# msg_deg
knockmsg <- read_rds('data/diff_expr_groups.rds')
msg_list <- read_lines('data/msg_list.txt')

# combine
diff_expr <- bind_rows(
    list(
        TF = select(knocktf, KD = TF, ID = Gene, FC = Log2FC, PVAL = P_value),
        MSG = select(knockmsg, KD = group, ID, FC = logFC, PVAL = P.Value)
    ),
    .id = 'Type')

diff_expr %>%
    filter(ID %in% unique(diff_expr$KD)) %>%
    mutate(relation2 = ifelse(FC > 0, -1, 1),
           Source = 'data',
           FC = round(FC, 2),
           PVAL = round(PVAL, 3)) %>%
    select(node1 = KD, 
           node2 = ID,
           relation2,
           everything()) %>%
    write_csv('data/data_interactions.csv')

# try merge
string_interactions <- read_csv('data/curated_interactions.csv')
data_interactions <- read_csv('data/data_interactions.csv')
dim(string_interactions) + dim(data_interactions)
full_join(string_interactions,
          data_interactions) %>%
    mutate(Source = paste(Source, source)) %>%
    write_csv('data/combined_interactions.csv')
