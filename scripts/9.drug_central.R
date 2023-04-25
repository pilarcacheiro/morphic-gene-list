### drugcentral information

# load libraries and auxiliary scripts ------------------------------------

if (!require("tidyverse")) install.packages("tidyverse")
library("tidyverse")

source("./scripts/auxiliary/hgnc_symbol_checker.R")


# import current gene list ------------------------------------------------

gene_list <- read_delim ("./data/processed/gene_list_reactome_annotations.txt") %>%
  dplyr::select(hgnc_id, gene_symbol, JAX, MSK, NWU, UCSF) %>%
  distinct()


# hgnc --------------------------------------------------------------------

hgnc_uniprot <- protein.coding.genes %>%
  select(hgnc_id, uniprot_ids)


# import drug central -----------------------------------------------------

drug_central <- read_delim ("https://unmtid-shinyapps.net/download/DrugCentral/2021_09_01/drug.target.interaction.tsv.gz") %>%
  filter(ORGANISM == "Homo sapiens") %>%
  select(GENE, TARGET_CLASS, TDL, DRUG_NAME) %>%
  distinct() %>%
  separate_rows(GENE,  sep = "\\|") %>%
  separate_rows(TDL,  sep = "\\|") %>%
  distinct() %>%
  group_by(GENE) %>%
  summarise(TARGET_CLASS = paste0(unique(TARGET_CLASS), collapse = "|"),
            TDL = paste0(unique(TDL), collapse = "|"),
            DRUG_NAME = paste0(unique(DRUG_NAME), collapse = "|")) %>%
  distinct()



# import TCRD -------------------------------------------------------------

tcrd <- read_delim ("http://juniper.health.unm.edu/tcrd/download/PharosTCRD_UniProt_Mapping.tsv") %>%
  select(UniProt_accession, TDL) %>%
  inner_join(hgnc_uniprot, by = c("UniProt_accession" = "uniprot_ids")) %>%
  select(hgnc_id, TDL) %>%
  distinct()


### symbols check

drug_symbols_check <-  hgnc.checker(drug_central$GENE, protein.coding.genes) %>%
  mutate(locus_type = "gene_with_protein_product") %>%
  filter(hgnc_id != "-") %>%
  select(hgnc_id, gene_symbol)

drug_central_hgnc_id <- drug_central %>%
  inner_join(drug_symbols_check, by = c("GENE" = "gene_symbol")) %>%
  select(-c(GENE, type)) %>%
  rename(drugcentral_target_class = TARGET_CLASS,
         drug_central_TDL = TDL,
         drug_central_drug_name = DRUG_NAME)


### gene list and drug central

gene_list_drug_central <- gene_list %>%
  left_join(tcrd) %>%
  left_join(drug_central_hgnc_id) %>%
  replace(is.na(.),"-") 


# export annotations ------------------------------------------------------

write.table(gene_list_drug_central,
            "./data/processed/drug_central_pharos.txt", quote = F,
            sep = "\t",row.names = F)