### panther families annotations

# load libraries and auxiliary scripts ------------------------------------

if (!require("tidyverse")) install.packages("tidyverse")
library("tidyverse")

if (!require("enrichplot")) install.packages("enrichplot")
library("enrichplot")

#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")

#BiocManager::install("ReactomePA")
library(ReactomePA)
library(reactome.db)

# import data -------------------------------------------------------------

gene_list <- read_delim ("./data/processed/gene_lists_merged.txt") %>%
  dplyr::select(hgnc_id, JAX, MSK, NWU, UCSF)


# retrieve only the annotations -------------------------------------------

genes_entrez <- read_delim("http://ftp.ebi.ac.uk/pub/databases/genenames/hgnc/tsv/locus_types/gene_with_protein_product.txt") %>%
  dplyr::select(hgnc_id, entrez_id)

gene_list_entrez <- gene_list %>%
  left_join(genes_entrez) %>%
  mutate(entrez_id = as.character(entrez_id))

path_name <- toTable(reactomePATHNAME2ID)

path_anot <- toTable(reactomeEXTID2PATHID) %>%
  filter(gene_id %in% gene_list_entrez$entrez_id) %>%
  dplyr::select(gene_id, DB_ID) %>%
  distinct() %>%
  inner_join(path_name) %>%
  dplyr::rename(reactome_path_id = DB_ID,
                reactome_path_name = path_name) %>%
  filter(grepl("Homo sapiens", reactome_path_name)) %>%
  mutate(reactome_path_name = stringr::str_sub(reactome_path_name, start = 15)) %>%
  distinct() %>%
  group_by(gene_id) %>%
  summarise(reactome_path_id = paste0(reactome_path_id, collapse = "|"),
            reactome_path_name = paste0(reactome_path_name, collapse = "|")) %>%
  replace(is.na(.),"-") 


gene_list_path <- gene_list_entrez %>%
  left_join(path_anot, by = c("entrez_id" = "gene_id"))%>%
  replace(is.na(.),"-")  %>%
  dplyr::select(-entrez_id)

hgnc_symbol <- read_delim("./data/processed/gene_lists_merged_impc_cells_data_disease.txt",
                          delim = "\t") %>%
  dplyr::select(hgnc_id, gene_symbol)

list_path_anot <- hgnc_symbol %>%
  left_join(gene_list_path) %>%
  relocate(hgnc_id, gene_symbol)


write.table(list_path_anot,
            "./data/processed/gene_list_reactome_annotations.txt", quote = F,
            sep = "\t",row.names = F)

# pathways enrichment -----------------------------------------------------

genes_universe <- read_delim("http://ftp.ebi.ac.uk/pub/databases/genenames/hgnc/tsv/locus_types/gene_with_protein_product.txt") %>%
  dplyr::select(hgnc_id, entrez_id)

jax <- gene_list %>%
  inner_join(genes_universe, by = c("hgnc_id" = "hgnc_id")) %>%
  filter(JAX == "JAX") %>%
  pull(entrez_id) %>%
  as.character(.)

msk <- gene_list %>%
  inner_join(genes_universe, by =  c("hgnc_id" = "hgnc_id")) %>%
  filter(MSK == "MSK") %>%
  pull(entrez_id) %>%
  as.character(.)

nwu <- gene_list %>%
  inner_join(genes_universe, by =  c("hgnc_id" = "hgnc_id")) %>%
  filter(NWU == "NWU") %>%
  pull(entrez_id) %>%
  as.character(.)

ucsf <- gene_list %>%
  inner_join(genes_universe, by =  c("hgnc_id" = "hgnc_id")) %>%
  filter(UCSF == "UCSF") %>%
  pull(entrez_id) %>%
  as.character(.)

universe <- as.character(genes_universe$entrez_id)

# jax ---------------------------------------------------------------------

jax_react <- enrichPathway(gene = jax, pvalueCutoff = 0.05, readable = TRUE)
head(jax_react)

jax_react_sim <- pairwise_termsim(jax_react)
emapplot(jax_react_sim)


# msk ---------------------------------------------------------------------

msk_react <- enrichPathway(gene = msk, pvalueCutoff = 0.05, readable = TRUE)
head(msk_react)

msk_react_sim <- pairwise_termsim(msk_react)
emapplot(msk_react_sim)

# nwu ---------------------------------------------------------------------

nwu_react <- enrichPathway(gene = nwu, pvalueCutoff = 0.05, readable = TRUE)
head(nwu_react)

nwu_react_sim <- pairwise_termsim(nwu_react)
emapplot(nwu_react_sim)


# ucsf ---------------------------------------------------------------------

ucsf_react <- enrichPathway(gene = ucsf, pvalueCutoff = 0.05, readable = TRUE)
head(ucsf_react)

ucsf_react_sim <- pairwise_termsim(ucsf_react)
emapplot(ucsf_react_sim)





