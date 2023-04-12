### merge gene lists
### check gene symbols

# load libraries and auxiliary scripts ------------------------------------

if (!require("tidyverse")) install.packages("tidyverse")
library("tidyverse")

source("./scripts/auxiliary/hgnc_symbol_checker.R")

# import gene lists, check symbols and merge ------------------------------

hgnc <- protein.coding.genes

jax <- read_delim("./data/raw/JAX gene list - Sheet1.tsv", delim = "\t") 

msk <- read_delim("./data/raw/MSK gene list - Sheet1.tsv", delim = "\t")

nwu <- read_delim("./data/raw/NWU gene list - Sheet1.tsv", delim = "\t")

ucsf <- read_delim("./data/raw/InitialGeneList_UCSF_MorPhiC.csv", delim = ",") %>%
  filter(!is.na(`GENE/Locus`))



### symbols check

jax_symbols_check <-  hgnc.checker(jax$`Gene Symbol`, hgnc) %>%
  mutate(locus_type = "gene_with_protein_product") %>%
  mutate(dpc = "JAX")

msk_symbols_check <-  hgnc.checker(msk$Gene_Symbol, hgnc) %>%
  filter(!is.na(gene_symbol)) %>%
  mutate(hgnc_id = ifelse(hgnc_id == "-", "HGNC:14575", hgnc_id)) %>%
  mutate(type = ifelse(type !="approved_symbol", "approved_symbol", type)) %>%
  mutate(locus_type = ifelse(hgnc_id == "HGNC:14575", "RNA_long_non_coding",
                             "gene_with_protein_product")) %>%
  mutate(dpc = "MSK")

nwu_symbols_check <-  hgnc.checker(nwu$Gene_Symbol, hgnc) %>%
  mutate(locus_type = "gene_with_protein_product") %>%
  mutate(dpc = "NWU")

ucsf_symbols_check <-  hgnc.checker(ucsf$`GENE/Locus`, hgnc) %>%
  mutate(hgnc_id = ifelse(hgnc_id == "-", "HGNC:28272", hgnc_id)) %>%
  mutate(type = ifelse(type =="Notfound.ProteinCoding.Symbol", "previous_symbol", type)) %>%
  mutate(locus_type = ifelse(hgnc_id == "HGNC:28272", "pseudogene",
                             "gene_with_protein_product")) %>%
  mutate(dpc = "UCSF")


### merge symbols

all_symbols_check <- jax_symbols_check %>%
  bind_rows(msk_symbols_check) %>%
  bind_rows(nwu_symbols_check) %>%
  bind_rows(ucsf_symbols_check)

all_symbols_check_wider <- all_symbols_check %>%
  pivot_wider(names_from = dpc, values_from = dpc) %>%
  replace(is.na(.),"-") %>%
  left_join(jax %>% dplyr::select(1,2), by = c("gene_symbol" = "Gene Symbol")) %>%
  dplyr::rename(JAX_MorPhiC_Code = "JAX MorPhiC Code") %>%
  left_join(msk, by = c("gene_symbol" = "Gene_Symbol")) %>%
  dplyr::rename(MSK_category = `Gene Category`) %>%
  relocate(gene_symbol) %>%
  replace(is.na(.),"-")


# export checked lists ----------------------------------------------------

write.table(all_symbols_check_wider, "./data/processed/gene_lists_merged.txt",
            quote = F, sep = "\t", row.names = F)
