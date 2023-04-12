### merge gene lists
### with komp/impc - mgi data

# load libraries and auxiliary scripts ------------------------------------

if (!require("tidyverse")) install.packages("tidyverse")
library("tidyverse")

source("./scripts/auxiliary/hgnc_symbol_checker.R")


# merge list, impc and mgi data -------------------------------------------

hgnc <- read_delim("http://ftp.ebi.ac.uk/pub/databases/genenames/hgnc/tsv/locus_types/gene_with_protein_product.txt") %>%
  dplyr::rename(gene_symbol = symbol, mgi_id = mgd_id) %>%
  dplyr::select(hgnc_id,mgi_id)%>%
  separate_rows(mgi_id,sep ="\\|") %>%
  filter(!is.na(mgi_id)) 

hgnc_dups <- unique(hgnc$hgnc_id[duplicated(hgnc$hgnc_id)])
mgi_dups <- unique(hgnc$mgi_id[duplicated(hgnc$mgi_id)])

hgnc_ortho_nodups <- hgnc %>%
  filter(!mgi_id %in% mgi_dups) %>%
  filter(!hgnc_id %in% hgnc_dups) %>%
  distinct()

hgnc_symbol <- read_delim("http://ftp.ebi.ac.uk/pub/databases/genenames/hgnc/tsv/locus_types/gene_with_protein_product.txt") %>% 
  dplyr::rename(gene_symbol = symbol) %>%
  dplyr::select(hgnc_id, gene_symbol) %>%
  filter(hgnc_id %in% hgnc_ortho_nodups$hgnc_id)

one2one <- hgnc_ortho_nodups

list <- read_delim("./data/processed/gene_lists_merged.txt", delim = "\t") 

impc_via <- read_delim("http://ftp.ebi.ac.uk/pub/databases/impc/all-data-releases/release-18.0/results/viability.csv.gz", 
                       delim = ",")  %>%
  dplyr::select('Gene Accession Id', 'Viability Phenotype HOMs/HEMIs', 'Comment') %>%
  dplyr::rename(mgi_id = 'Gene Accession Id',
         impc_viability = 'Viability Phenotype HOMs/HEMIs',
         impc_viability_notes = 'Comment') %>%
  distinct() %>%
  arrange(mgi_id, impc_viability) %>%
  group_by(mgi_id) %>%
  summarise(impc_viability = paste0(unique(impc_viability), collapse = "|"),
            impc_viability_notes = paste0(unique(impc_viability_notes), collapse = "|")) %>%
replace(is.na(.),"-") %>%
  mutate(impc_viability_notes = ifelse(impc_viability_notes == "NA","-",
                                       impc_viability_notes))

impc_pheno  <- read_delim("http://ftp.ebi.ac.uk/pub/databases/impc/all-data-releases/release-18.0/results/genotype-phenotype-assertions-ALL.csv.gz", 
                          delim = ",") %>%
  dplyr::select(marker_accession_id, zygosity,mp_term_name) %>%
  distinct() %>%
  dplyr::rename(mgi_id = marker_accession_id, 
         impc_zygosity = zygosity, 
         impc_phenotype = mp_term_name) %>%
  group_by(mgi_id, impc_zygosity) %>%
  summarise(impc_phenotypes = paste0(unique(impc_phenotype), collapse = "|")) %>%
  pivot_wider(names_from = impc_zygosity, values_from = impc_phenotypes) %>%
  replace(is.na(.),"-") %>%
  dplyr::rename(impc_phenotypes_homozygote = homozygote,
         impc_phenotypes_heterozygote  = heterozygote,
         impc_phenotypes_hemizygote = hemizygote) %>%
  replace(is.na(.),"-")

### merge with lists

list_impc <- list %>%
  left_join(one2one) %>%
  left_join(impc_via) %>%
  left_join(impc_pheno)  %>%
  replace(is.na(.),"-") %>%
  mutate(one2one_ortholog_hgnc = ifelse(mgi_id == "-", "-","y")) %>%
  dplyr::select(gene_symbol:MSK_category, one2one_ortholog_hgnc, mgi_id,
         impc_viability, impc_phenotypes_homozygote:impc_phenotypes_hemizygote) %>%
  replace(is.na(.),"-") 


# export checked lists ----------------------------------------------------

write.table(list_impc, "./data/processed/gene_lists_merged_impc_data.txt",
            quote = F, sep = "\t", row.names = F)