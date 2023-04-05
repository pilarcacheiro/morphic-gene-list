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

# import data -------------------------------------------------------------

gene_list <- read_delim ("./data/processed/gene_lists_merged.txt") %>%
  dplyr::select(gene_symbol, JAX, MSK, NWU, UCSF)

genes_universe <- read_delim("http://ftp.ebi.ac.uk/pub/databases/genenames/hgnc/tsv/locus_types/gene_with_protein_product.txt") %>%
  dplyr::select(symbol, entrez_id)

jax <- gene_list %>%
  inner_join(genes_universe, by = c("gene_symbol" = "symbol")) %>%
  filter(JAX == "JAX") %>%
  pull(entrez_id) %>%
  as.character(.)

msk <- gene_list %>%
  inner_join(genes_universe, by = c("gene_symbol" = "symbol")) %>%
  filter(MSK == "MSK") %>%
  pull(entrez_id) %>%
  as.character(.)

nwu <- gene_list %>%
  inner_join(genes_universe, by = c("gene_symbol" = "symbol")) %>%
  filter(NWU == "NWU") %>%
  pull(entrez_id) %>%
  as.character(.)

ucsf <- gene_list %>%
  inner_join(genes_universe, by = c("gene_symbol" = "symbol")) %>%
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





