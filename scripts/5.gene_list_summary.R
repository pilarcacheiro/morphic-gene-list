### descriptive analysis

# load libraries and auxiliary scripts ------------------------------------

if (!require("tidyverse")) install.packages("tidyverse")
library("tidyverse")

if (!require("UpSetR")) install.packages("UpSetR")
library("UpSetR")

if (!require("cowplot")) install.packages("cowplot")
library("cowplot")

source("./scripts/auxiliary/hgnc_symbol_checker.R")


# import merged list

gene_list <- read_delim("./data/processed/gene_lists_merged_impc_cells_data_disease.txt",
                        delim = "\t")

# overlap UpsetR

gene_upset <- gene_list %>%
  dplyr::select(gene_symbol,JAX, MSK,NWU, UCSF) %>%
  mutate(JAX  = ifelse(JAX == "JAX",1,0))%>%
  mutate(MSK  = ifelse(MSK == "MSK",1,0))%>%
  mutate(NWU  = ifelse(NWU == "NWU",1,0)) %>%
  mutate(UCSF  = ifelse(UCSF == "UCSF",1,0))

gene_upset <- as.data.frame(gene_upset)

upset(gene_upset, nsets = 4, keep.order = TRUE, point.size = 4, line.size = 1.5,
      mainbar.y.label = "Number of genes", sets.x.label = "Genes per centre", 
      text.scale = c(2, 2, 2, 2, 2, 2),
      set_size.show = TRUE)



# mouse lethal lines ------------------------------------------------------

## genes with orthologues

gene_list_mouse <- gene_list %>%
  filter(mgi_id !="-") %>%
  filter(impc_viability !="-") %>%
  mutate(impc_viability_2 = ifelse(!impc_viability %in% c("lethal","subviable","viable"),
                                   "conflicting", impc_viability)) %>%
  group_by(impc_viability_2) %>%
  tally() %>%
  mutate(impc_viability_3 = factor(impc_viability_2,
                                  levels = c("lethal","subviable","viable","conflicting"))) %>%
  mutate(percentage = (n/sum(n)*100))


plot_via_all <- ggplot(gene_list_mouse, aes(x=impc_viability_3, y=percentage)) +
  geom_bar(stat="identity", fill = "#E9D8A6") +
  theme_minimal() +
  ylab("% of genes") +
  xlab("IMPC preweaning viability assessment") +
  theme(legend.position="none") +
  theme(axis.title.y = 
          element_text(size = 10,margin = margin(t = 0, r = 15, b = 0, l = 0))) +
  theme(axis.title.x = 
          element_text(size = 10 ,margin = margin(t = 15, r = 0, b = 0, l = 0))) + 
  theme(axis.text.x = element_text(size = 11, face = "bold")) +
  theme(axis.text.y = element_text(size = 10)) +
  geom_text(data = gene_list_mouse,
            aes(label=paste0(round(percentage), " %")), 
            size = 3.5,nudge_y = 3)


## by center


gene_mouse <- gene_list %>%
  filter(mgi_id !="-") %>%
  filter(impc_viability !="-") %>%
  mutate(impc_viability_2 = ifelse(!impc_viability %in% c("lethal","subviable","viable"),
                                   "conflicting", impc_viability))  %>%
  dplyr::select(hgnc_id, impc_viability_2 )


gene_list_mouse_long  <- gene_list %>%
  dplyr::select(hgnc_id, JAX, MSK, NWU, UCSF) %>%
  pivot_longer(!hgnc_id, names_to = "center", values_to = "selected") %>%
  filter(selected !="-") %>%
  dplyr::select(hgnc_id, selected) %>%
  inner_join(gene_mouse)  %>%
  distinct() %>%
  dplyr::rename(center = selected)


gene_list_mouse_long_center <- gene_list_mouse_long %>%
  group_by(center) %>%
  tally() %>%
  dplyr::rename(n_center = n)

gene_list_mouse_long_center_via <- gene_list_mouse_long %>%
  group_by(center, impc_viability_2) %>%
  tally() %>%
  inner_join(gene_list_mouse_long_center) %>%
  mutate(percentage = (n / n_center)*100) %>%
  mutate(impc_viability_3 = factor(impc_viability_2,
                                   levels = c("lethal","subviable","viable","conflicting")))

plot_via_center <- ggplot(gene_list_mouse_long_center_via, 
       aes(x = impc_viability_3, y = percentage, fill = center)) +
  geom_bar(stat="identity") +
  facet_wrap(~center, nrow = 2) +
  theme_minimal() +
  ylab("% of genes") +
  xlab("IMPC preweaning viability assessment") +
  theme(legend.position="none") +
  theme(axis.title.y = 
          element_text(size = 10,margin = margin(t = 0, r = 15, b = 0, l = 0))) +
  theme(axis.title.x = 
          element_text(size = 10 ,margin = margin(t = 15, r = 0, b = 0, l = 0))) + 
  theme(axis.text.x = element_text(size = 11, face = "bold")) +
  theme(axis.text.y = element_text(size = 10)) +
  geom_text(data = gene_list_mouse_long_center_via,
            aes(label=paste0(round(percentage), " %")), 
            size = 3.5,nudge_y = 8) +
  scale_fill_manual(values=c("#0A9396","#EE9B00","#9B2226","#94D2BD"))

  
plot_grid(plot_via_all, plot_via_center, nrow = 1)


# cell line essentiality --------------------------------------------------

## genes with cell data

gene_list_cells <- gene_list %>%
  dplyr::select(hgnc_id, JAX, MSK, NWU, UCSF, cellular_essential) %>%
  filter(cellular_essential !="-") %>%
  mutate(cellular_essential = recode(cellular_essential, y = "Essential", n = "Non-essential"))

gene_score <- gene_list_cells %>%
  dplyr::select(hgnc_id, cellular_essential)

gene_list_cells_long  <- gene_list_cells %>%
  dplyr::select(hgnc_id, JAX, MSK, NWU, UCSF) %>%
  pivot_longer(!hgnc_id, names_to = "center", values_to = "selected") %>%
  filter(selected !="-") %>%
  dplyr::select(hgnc_id, selected) %>%
  inner_join(gene_score)  %>%
  distinct() %>%
  dplyr::rename(center = selected)

gene_score_count <- gene_score %>%
  group_by(cellular_essential) %>%
  tally() %>%
  mutate(percentage = (n/sum(n)*100))


gene_list_score_long_all <- gene_list_cells_long  %>%
  group_by(center) %>%
  tally() %>%
  dplyr::rename(n_center = n)

gene_list_score_long_center_count <- gene_list_cells_long %>%
  group_by(center, cellular_essential) %>%
  tally() %>%
  inner_join(gene_list_score_long_all) %>%
  mutate(percentage = (n / n_center)*100) 

plot_cells_all <- ggplot(gene_score_count, aes(x = cellular_essential, y = percentage)) +
  geom_bar(stat="identity", fill = "#E9D8A6") +
  theme_minimal() +
  ylab("% of genes") +
  xlab("Cellular essentiality") +
  theme(legend.position="none") +
  theme(axis.title.y = 
          element_text(size = 10,margin = margin(t = 0, r = 15, b = 0, l = 0))) +
  theme(axis.title.x = 
          element_text(size = 10 ,margin = margin(t = 15, r = 0, b = 0, l = 0))) + 
  theme(axis.text.x = element_text(size = 11, face = "bold")) +
  theme(axis.text.y = element_text(size = 10)) +
  geom_text(data = gene_score_count,
            aes(label=paste0(round(percentage), " %")), 
            size = 3.5,nudge_y = 3)
    
    
plot_cells_center <- ggplot(gene_list_score_long_center_count, 
                          aes(x = cellular_essential, y = percentage, fill = center)) +
  geom_bar(stat="identity") +
  facet_wrap(~center, nrow = 2) +
  theme_minimal() +
  ylab("% of genes") +
  xlab("Cellular essentiality") +
  theme(legend.position="none") +
  theme(axis.title.y = 
          element_text(size = 10,margin = margin(t = 0, r = 15, b = 0, l = 0))) +
  theme(axis.title.x = 
          element_text(size = 10 ,margin = margin(t = 15, r = 0, b = 0, l = 0))) + 
  theme(axis.text.x = element_text(size = 11, face = "bold")) +
  theme(axis.text.y = element_text(size = 10)) +
  geom_text(data = gene_list_score_long_center_count,
            aes(label=paste0(round(percentage), " %")), 
            size = 3.5,nudge_y = 8) +
  scale_fill_manual(values=c("#0A9396","#EE9B00","#9B2226","#94D2BD"))


plot_grid(plot_cells_all, plot_cells_center, nrow = 1)
    
    

# disease genes -----------------------------------------------------------

## genes associated to disease

gene_list_disease <- gene_list %>%
  dplyr::select(hgnc_id, JAX, MSK, NWU, UCSF, omim, dd) %>%
  mutate(omim_gene = recode(omim, "-" = "nonOMIM", omim_disease = "OMIM")) %>%
  mutate(omim_gene = factor(omim_gene, levels = c("OMIM","nonOMIM")))



gene_omim <- gene_list_disease%>%
  dplyr::select(hgnc_id, omim_gene)

gene_list_omim_long  <- gene_list_disease %>%
  dplyr::select(hgnc_id, JAX, MSK, NWU, UCSF) %>%
  pivot_longer(!hgnc_id, names_to = "center", values_to = "selected") %>%
  filter(selected !="-") %>%
  dplyr::select(hgnc_id, selected) %>%
  inner_join(gene_omim )  %>%
  distinct() %>%
  dplyr::rename(center = selected)

gene_omim_count <- gene_omim%>%
  group_by(omim_gene) %>%
  tally() %>%
  mutate(percentage = (n/sum(n)*100))


gene_list_omim_long_all <- gene_list_omim_long   %>%
  group_by(center) %>%
  tally() %>%
  dplyr::rename(n_center = n)

gene_list_omim_long_center_count <- gene_list_omim_long %>%
  group_by(center, omim_gene) %>%
  tally() %>%
  inner_join(gene_list_omim_long_all) %>%
  mutate(percentage = (n / n_center)*100) 

plot_omim_all <- ggplot(gene_omim_count, aes(x = omim_gene, y = percentage)) +
  geom_bar(stat="identity", fill = "#E9D8A6") +
  theme_minimal() +
  ylab("% of genes") +
  xlab("Mendelian disease association") +
  theme(legend.position="none") +
  theme(axis.title.y = 
          element_text(size = 10,margin = margin(t = 0, r = 15, b = 0, l = 0))) +
  theme(axis.title.x = 
          element_text(size = 10 ,margin = margin(t = 15, r = 0, b = 0, l = 0))) + 
  theme(axis.text.x = element_text(size = 11, face = "bold")) +
  theme(axis.text.y = element_text(size = 10)) +
  geom_text(data = gene_omim_count,
            aes(label=paste0(round(percentage), " %")), 
            size = 3.5,nudge_y = 3)


plot_omim_center <- ggplot(gene_list_omim_long_center_count, 
                            aes(x = omim_gene, y = percentage, fill = center)) +
  geom_bar(stat="identity") +
  facet_wrap(~center, nrow = 2) +
  theme_minimal() +
  ylab("% of genes") +
  xlab("Mendelian disease association") +
  theme(legend.position="none") +
  theme(axis.title.y = 
          element_text(size = 10,margin = margin(t = 0, r = 15, b = 0, l = 0))) +
  theme(axis.title.x = 
          element_text(size = 10 ,margin = margin(t = 15, r = 0, b = 0, l = 0))) + 
  theme(axis.text.x = element_text(size = 11, face = "bold")) +
  theme(axis.text.y = element_text(size = 10)) +
  geom_text(data = gene_list_omim_long_center_count,
            aes(label=paste0(round(percentage), " %")), 
            size = 3.5,nudge_y = 8) +
  scale_fill_manual(values=c("#0A9396","#EE9B00","#9B2226","#94D2BD"))


plot_grid(plot_omim_all, plot_omim_center, nrow = 1)


 
    
