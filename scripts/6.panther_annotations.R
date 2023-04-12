### panther families annotations


### we need to bring the panther ontology to be able to compare
### protein familiies at the same level that are retrieved through
### the bioconductor package
### http://data.pantherdb.org/PANTHER17.0/ontology/Protein_Class_17.0
### http://data.pantherdb.org/PANTHER17.0/ontology/Protein_class_relationship

### currently the famlies are retrieved through the browser using the uniprotids
### from gene_uniprotid.txt file

# load libraries and auxiliary scripts ------------------------------------

if (!require("tidyverse")) install.packages("tidyverse")
library("tidyverse")


#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")

#BiocManager::install("PANTHER.db")
library("PANTHER.db")


# import data -------------------------------------------------------------


human_uniprot <-  read_delim("http://ftp.ebi.ac.uk/pub/databases/genenames/hgnc/tsv/locus_types/gene_with_protein_product.txt")  %>%
  dplyr::select(hgnc_id, uniprot_ids) %>%
  dplyr::rename(uniprot_id = uniprot_ids) %>%
  dplyr::distinct()

write.table(human_uniprot, "./data/raw/gene_uniprotid.txt",
quote = F, sep = "\t",row.names = F)

gene_list <- read_delim("./data/processed/gene_lists_merged_impc_cells_data_disease.txt",
                        delim = "\t")

panther <- read_delim("./data/processed/pantherGeneList.txt",
                      delim = "\t",  col_names = FALSE) %>%
  dplyr::select("X2","X5") %>%
  dplyr::rename(uniprot_id = X2, protein_class = X5)

hgnc_symbol <- read_delim("./data/processed/gene_lists_merged_impc_cells_data_disease.txt",
                          delim = "\t") %>%
  dplyr::select(hgnc_id, gene_symbol)


# extract info from panther -----------------------------------------------

columns(PANTHER.db)
keytypes(PANTHER.db)
head(keys(PANTHER.db))

pthOrganisms(PANTHER.db) <- "HUMAN"

cols <- c("UNIPROT","CLASS_ID","CLASS_TERM")

cols <- c("UNIPROT","FAMILY_ID","FAMILY_TERM")

keys <- human_uniprot$uniprot_id

kt <- "UNIPROT"

#select(PANTHER.db, keys, cols, kt)

#panther_annotation <- select(PANTHER.db, keys = "Q9UNA3", columns = cols, keytype = kt)

#panther_class <- panther_annotation %>%
#  dplyr::select(1,2,3) %>%
#  dplyr::rename(uniprot_id = UNIPROT,
#                class_id = CLASS_ID,
#                class_term = CLASS_TERM) %>%
#  dplyr::inner_join(human_uniprot) %>%
#  dplyr::relocate(hgnc_id) %>%
#  dplyr::distinct()

term <- "PC00220"
select(PANTHER.db,term, "CLASS_TERM","CLASS_ID")
ancestors <- traverseClassTree(PANTHER.db,term,scope="ANCESTOR")

### we need to bring the panther ontology to be able to compare
### protein familiies at the same level
### http://data.pantherdb.org/PANTHER17.0/ontology/Protein_Class_17.0
### http://data.pantherdb.org/PANTHER17.0/ontology/Protein_class_relationship


# file from browser -------------------------------------------------------


gene_list_panther <- gene_list %>%
  inner_join(human_uniprot) %>%
  inner_join(panther ) %>%
  filter(!is.na(protein_class))%>%
  dplyr::select(hgnc_id, JAX, MSK, NWU, UCSF, protein_class) 



# annotations to export individually --------------------------------------


gene_list_panther_hgnc <- hgnc_symbol %>%
  left_join(gene_list_panther) %>%
  relocate(hgnc_id, gene_symbol)


write.table(gene_list_panther_hgnc,
            "./data/processed/gene_list_panther_protein_class.txt", quote = F,
            sep = "\t",row.names = F)

gene_list_panther_export <- gene_list %>%
  inner_join(human_uniprot) %>%
  inner_join(panther ) %>%
  filter(!is.na(protein_class))%>%
  dplyr::select(uniprot_id, JAX, MSK, NWU, UCSF)

write.table(gene_list_panther_export,
            "./data/processed/uniprot_centers.txt", quote = F,
            sep = "\t",row.names = F)


# protein class analysis --------------------------------------------------



gene_panther  <- gene_list_panther %>%
  dplyr::select(hgnc_id, protein_class)

gene_list_panther_long  <- gene_list_panther %>%
  dplyr::select(hgnc_id, JAX, MSK, NWU, UCSF) %>%
  pivot_longer(!hgnc_id, names_to = "center", values_to = "selected") %>%
  filter(selected !="-") %>%
  dplyr::select(hgnc_id, selected) %>%
  inner_join(gene_panther)  %>%
  distinct() %>%
  dplyr::rename(center = selected)

gene_panther_count <- gene_panther %>%
  group_by(protein_class) %>%
  tally() %>%
  mutate(percentage = (n/805 *100)) %>%
  arrange(-percentage)


gene_list_panther_long_all <- gene_list_panther_long %>%
  group_by(center) %>%
  tally() %>%
  dplyr::rename(n_center = n)

gene_list_panther_long_center_count <- gene_list_panther_long %>%
  group_by(center, protein_class) %>%
  tally() %>%
  inner_join(gene_list_panther_long_all) %>%
  mutate(percentage = (n / n_center)*100) 

gene_panther_count_b <- gene_panther_count %>%
  dplyr::slice(1:10) %>%
  mutate(protein_class = factor(protein_class, levels = c(gene_panther_count$protein_class[10:1])))

gene_list_panther_long_center_count_b <- gene_list_panther_long_center_count %>%
  arrange(desc(percentage)) %>% 
  group_by(center) %>% 
  dplyr::slice(1:10)

plot_panther_all <- ggplot(gene_panther_count_b, aes(x = protein_class, y = percentage)) +
  geom_bar(stat="identity", fill = "#E9D8A6") +
  coord_flip() +
  theme_minimal() +
  ylab("% of genes") +
  xlab("Protein class") +
  theme(legend.position="none") +
  theme(axis.title.y = 
          element_text(size = 10,margin = margin(t = 0, r = 15, b = 0, l = 0))) +
  theme(axis.title.x = 
          element_text(size = 10 ,margin = margin(t = 15, r = 0, b = 0, l = 0))) + 
  theme(axis.text.x = element_text(size = 10, )) +
  theme(axis.text.y = element_text(size = 10)) +
  geom_text(data = gene_panther_count_b,
            aes(label=paste0(round(percentage), " %")), 
            size = 3.5,nudge_y = 3)


plot_panther_center <- ggplot(gene_list_panther_long_center_count_b, 
                           aes(x = protein_class, y = percentage, fill = center)) +
  geom_bar(stat="identity") +
  coord_flip() +
  facet_wrap(~center, nrow = 1) +
  theme_minimal() +
  ylab("% of genes") +
  xlab("Protein class") +
  theme(legend.position="none") +
  theme(axis.title.y = 
          element_text(size = 10,margin = margin(t = 0, r = 15, b = 0, l = 0))) +
  theme(axis.title.x = 
          element_text(size = 10 ,margin = margin(t = 15, r = 0, b = 0, l = 0))) + 
  theme(axis.text.x = element_text(size = 11)) +
  theme(axis.text.y = element_text(size = 10)) +
  geom_text(data = gene_list_panther_long_center_count_b,
            aes(label=paste0(round(percentage), " %")), 
            size = 3.5,nudge_y = 8) +
  scale_fill_manual(values=c("#0A9396","#EE9B00","#9B2226","#94D2BD"))


plot_grid(plot_panther_all, plot_panther_center, nrow = 1, rel_widths = c(0.5,1))




