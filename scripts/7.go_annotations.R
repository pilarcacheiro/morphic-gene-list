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




# retrieve only the annotations -------------------------------------------


genes_entrez <- read_delim("http://ftp.ebi.ac.uk/pub/databases/genenames/hgnc/tsv/locus_types/gene_with_protein_product.txt") %>%
  dplyr::select(hgnc_id, entrez_id) %>%
  mutate(entrez_id = as.character(entrez_id))

gene_list_entrez <- gene_list %>%
  left_join(genes_entrez)

go_anot <- toTable(org.Hs.egGO) %>%
  filter(gene_id %in% gene_list_entrez$entrez_id) %>%
  dplyr::select(gene_id, go_id, Ontology) %>%
  distinct()

keys <- unique(go_anot$go_id)

go_terms <- AnnotationDbi::select(GO.db, keys = keys, keytype="GOID", columns=c("TERM") ) %>%
  dplyr::rename(go_id = GOID,
                go_term = TERM)

go_anot_term <- go_anot %>%
  left_join(go_terms) %>%
  group_by(gene_id, Ontology) %>%
  summarise(go_ids = paste0(go_id, collapse = "|"),
            go_terms = paste0(go_term, collapse = "|")) %>%
  pivot_wider(names_from = Ontology,
              values_from = c("go_ids","go_terms")) %>%
  dplyr::select(gene_id, go_ids_BP, go_terms_BP,
                go_ids_MF, go_terms_MF,
                go_ids_CC, go_terms_CC)%>%
  replace(is.na(.),"-")


hgnc_symbol <- read_delim("./data/processed/gene_lists_merged_impc_cells_data_disease.txt",
                          delim = "\t") %>%
  dplyr::select(hgnc_id, gene_symbol)

list_go_anot <- gene_list_entrez %>%
  left_join(go_anot_term, by = c("entrez_id" = "gene_id")) %>%
  dplyr::select(-entrez_id) %>%
  replace(is.na(.),"-")

gene_list_go_hgnc <- hgnc_symbol %>%
  left_join(list_go_anot) %>%
  relocate(hgnc_id, gene_symbol)

write.table(gene_list_go_hgnc,
            "./data/processed/gene_list_go_annotations.txt", quote = F,
            sep = "\t",row.names = F)


# enrichment analysis -----------------------------------------------------



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


