### go annotations and enrichment

# load libraries and auxiliary scripts ------------------------------------

if (!require("tidyverse")) install.packages("tidyverse")
library("tidyverse")

#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")

#BiocManager::install("clusterProfiler")
library("clusterProfiler")


#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")

#BiocManager::install("org.Hs.eg.db")
library("org.Hs.eg.db")

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

jax_go <- enrichGO(gene          = jax,
               universe      = universe,
               OrgDb         = org.Hs.eg.db,
               ont           = "BP",
               pAdjustMethod = "BH",
               pvalueCutoff  = 0.01,
               qvalueCutoff  = 0.05,
               readable      = TRUE)

jax_go_df <- jax_go@result %>%
  filter(Count > 5) %>%
  filter(p.adjust < 0.01)

jax_go_revigo  <- jax_go_df %>%
  dplyr::select(ID,p.adjust) %>%
  dplyr::rename(go_id = ID, p_adjust = p.adjust)


# msk ---------------------------------------------------------------------

msk_go <- enrichGO(gene          = msk,
                   universe      = universe,
                   OrgDb         = org.Hs.eg.db,
                   ont           = "BP",
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.01,
                   qvalueCutoff  = 0.05,
                   readable      = TRUE)

msk_go_df <- msk_go@result %>%
  filter(Count > 5) %>%
  filter(p.adjust < 0.01)

msk_go_revigo  <- msk_go_df %>%
  dplyr::select(ID,p.adjust) %>%
  dplyr::rename(go_id = ID, p_adjust = p.adjust)


# nwu ---------------------------------------------------------------------

nwu_go <- enrichGO(gene          = nwu,
                   universe      = universe,
                   OrgDb         = org.Hs.eg.db,
                   ont           = "BP",
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.01,
                   qvalueCutoff  = 0.05,
                   readable      = TRUE)

nwu_go_df <- nwu_go@result %>%
  filter(Count > 5) %>%
  filter(p.adjust < 0.01)

nwu_go_revigo  <- nwu_go_df %>%
  dplyr::select(ID,p.adjust) %>%
  dplyr::rename(go_id = ID, p_adjust = p.adjust)


# ucsf ---------------------------------------------------------------------

ucsf_go <- enrichGO(gene         = ucsf,
                   universe      = universe,
                   OrgDb         = org.Hs.eg.db,
                   ont           = "BP",
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.01,
                   qvalueCutoff  = 0.05,
                   readable      = TRUE)

ucsf_go_df <- ucsf_go@result %>%
  filter(Count > 5) %>%
  filter(p.adjust < 0.01)

ucsf_go_revigo  <- ucsf_go_df %>%
  dplyr::select(ID, p.adjust) %>%
  dplyr::rename(go_id = ID, p_adjust = p.adjust)


# export files for REVIGO -------------------------------------------------

write.table(jax_go_revigo, "./data/processed/jax_go_pvalues.txt",
            quote = F, sep = "\t", row.names = F)

write.table(msk_go_revigo, "./data/processed/msk_go_pvalues.txt",
            quote = F, sep = "\t", row.names = F)

write.table(nwu_go_revigo, "./data/processed/nwu_go_pvalues.txt",
            quote = F, sep = "\t", row.names = F)

write.table(ucsf_go_revigo, "./data/processed/ucsf_go_pvalues.txt",
            quote = F, sep = "\t", row.names = F)


