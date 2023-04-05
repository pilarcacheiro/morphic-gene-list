### merge gene lists
### with disease annotations

# load libraries and auxiliary scripts ------------------------------------

if (!require("tidyverse")) install.packages("tidyverse")
library("tidyverse")

source("./scripts/auxiliary/hgnc_symbol_checker.R")


# cell lines --------------------------------------------------------------

list_impc <- read_delim("./data/processed/gene_lists_merged_impc_cells_data.txt",
                        delim = "\t")

# omim genes --------------------------------------------------------------

omim <- read_delim("./data/raw/morbidmap_complete.txt", delim = "\t") %>%
  select(hgnc_id, omim_title) %>%
  distinct() %>%
  filter(hgnc_id !="-") %>%
  filter(!is.na(hgnc_id)) %>%
  filter(!is.na(omim_title)) %>%
  group_by(hgnc_id) %>%
  summarise(omim = "omim_disease",
            omim_name = paste0(omim_title, collapse = "|"))


# DD genes ----------------------------------------------------------------

dd <- read_delim("./data/raw/DDG2P_7_3_2023.csv.gz", delim = ",") %>%
  rename(hgnc = "hgnc id") %>%
  mutate(hgnc_id = paste0("HGNC:", hgnc)) %>%
  mutate(dd = "dd" ) %>%
  rename(dd_name = "disease name",
         dd_confidence_category = "confidence category",
         dd_allelic_requirement = "allelic requirement") %>%
  select(hgnc_id, dd, dd_name, dd_confidence_category, dd_allelic_requirement) %>%
  group_by(hgnc_id) %>%
  summarise(dd = unique(dd), 
            dd_name = paste0(dd_name, collapse = "|"),
            dd_confidence_category = paste0(dd_confidence_category, collapse = "|"),
            dd_allelic_requirement = paste0(dd_allelic_requirement, collapse = "|"))

# merge mouse cell and disease --------------------------------------------


mouse_cells_disease <- mouse_cells %>%
  left_join(omim) %>%
  mutate(omim = ifelse(is.na(omim),"-", omim)) %>%
  mutate(omim_name = ifelse(is.na(omim_name),"-", omim_name)) %>%
  left_join(dd) %>%
  mutate(dd = ifelse(is.na(dd),"-", dd)) %>%
  mutate(dd_name = ifelse(is.na(dd_name),"-", dd_name)) %>%
  mutate(dd_confidence_category = ifelse(is.na(dd_confidence_category),"-", dd_confidence_category)) %>%
  mutate(dd_allelic_requirement = ifelse(is.na(dd_allelic_requirement),"-", dd_allelic_requirement))


# export current data -----------------------------------------------------

write.table(mouse_cells_disease, "./data/processed/gene_lists_merged_impc_cells_data_disease.txt",
            quote = F, sep = "\t", row.names = F)



