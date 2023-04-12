### merge gene lists
### check gene symbols
### currently depmap file too big for github (even gz)
### downoad from https://depmap.org/portal/download/all/
### 22Q4
### CRISPRGeneEffect.csv

# load libraries and auxiliary scripts ------------------------------------

if (!require("tidyverse")) install.packages("tidyverse")
library("tidyverse")

source("./scripts/auxiliary/hgnc_symbol_checker.R")


# cell lines --------------------------------------------------------------

list_impc <- read_delim("./data/processed/gene_lists_merged_impc_data.txt",
                        delim = "\t")

### Model.csv file from DepMap

models <- read_delim("./data/raw/Model.csv", delim = ",") 

unique(models$OncotreeLineage)

cns_brain <- models %>% filter(OncotreeLineage == "CNS/Brain")


cell0 <- read_delim("./data/raw/CRISPRGeneEffect.csv", delim = ",") 
colnames(cell0) <-  str_split(names(cell0), "\\ ", simplify=T)[,1]

colnames(cell0)[colnames(cell0) == "...1"] <- "model_id"


cell1  = cell0[,-1]

cell2 <- t(cell1)
cell3  <- as.data.frame(cell2)
colnames(cell3) = cell0$model_id

cell4 = cell3

cell4[cell4 >= -0.5] <- 0

cell4[cell4 != 0] <- 1

cell5 <- cell4 %>%
  mutate(n_essential = rowSums(.)) %>%
  mutate(percentage_essential = (n_essential/1078)*100)

cell6 <- cell4 %>% 
  dplyr::select(matches(cns_brain$ModelID))%>%
  mutate(n_essential_cns_brain = rowSums(.)) %>%
  mutate(percentage_essential_cns_brain = (n_essential_cns_brain/93)*100)


cell7 <- cell3 %>%
  mutate(mean_score_all = rowMeans(.)) %>%
  mutate(gene_symbol = row.names(cell3)) %>%
  dplyr::select(gene_symbol, mean_score_all) #%>%
#  `rownames <-`( NULL )
rownames(cell7) <- NULL


cell8 <- cell3 %>% 
  dplyr::select(matches(cns_brain$ModelID))%>%
  mutate(mean_score_cns_brain = rowMeans(.)) %>%
  mutate(gene_symbol = row.names(cell3)) %>%
  dplyr::select(gene_symbol, mean_score_cns_brain) #%>%
  #`rownames <-`( NULL )
rownames(cell8) <- NULL



cells_n <- cell5 %>%
  mutate(gene_symbol = row.names(cell5)) %>%
  dplyr::select(gene_symbol, n_essential, percentage_essential) %>%
  inner_join(cell6 %>%
               mutate(gene_symbol = row.names(cell6)) %>%
               dplyr::select(gene_symbol, n_essential_cns_brain, percentage_essential_cns_brain)) %>%
  dplyr::rename(n_essential_all_1078 = n_essential,
         percentage_essential_all = percentage_essential,
         n_essential_cns_brain_93 = n_essential_cns_brain) %>%
  inner_join(cell7) %>%
  inner_join(cell8) 


cells_depmap <- cells_n %>%
  dplyr::select(gene_symbol, percentage_essential_all,mean_score_all) %>%
  mutate(cell_essential =cut(mean_score_all, breaks=c(-Inf, -0.5, Inf), labels=c("y","n")))


# merge mouse and cells ---------------------------------------------------

mouse_cells <- list_impc %>%
  left_join(cells_depmap, by = c("gene_symbol" = "gene_symbol")) %>%
  mutate(cell_essential = as.character(cell_essential)) %>%
  mutate(cellular_essential = replace_na(cell_essential, "-")) %>%
  dplyr::select(-cell_essential)




# export checked lists with mouse and cell data ---------------------------

write.table(mouse_cells , "./data/processed/gene_lists_merged_impc_cells_data.txt",
            quote = F, sep = "\t", row.names = F)
