################################################################################
## Script: hgnc_symbol_checker.R
################################################################################

# required packages -------------------------------------------------------

if (!require("magrittr")) install.packages("magrittr")
library("magrittr")

if (!require("dplyr")) install.packages("dplyr")
library("dplyr")

if (!require("tidyr")) install.packages("tidyr")
library("tidyr")


# hgnc gene file ----------------------------------------------------------

# If the file doesn't already exist, we will read it via FTP
filename <- "gene_with_protein_product.txt"
if (!file.exists(filename)) {
  filename <- paste("ftp://ftp.ebi.ac.uk/pub/databases/genenames/hgnc/tsv/locus_types",
                    filename,
                    sep = "/")
}

protein.coding.genes <- readr::read_delim(filename,
                                          delim = "\t",
                                          col_names = TRUE)

## a given gene symbol can belong to one of different categories according to
## hgnc:
##   - symbol: The official gene symbol approved by the HGNC, which is
##             typically a short form of the gene name
##   - alias symbol: Alternative names for the gene. Aliases may be from
##                   literature, from other databases or may be added to
##                   represent membership of a gene group
##   - previous symbol: This field displays any symbols that were previously
##                      HGNC-approved nomenclature



# hgnc checker function ---------------------------------------------------


hgnc.checker <- function(gene_symbols, gene.file) {
  ## Purpose:
  ##   check gene symbols and retrieve hgnc id (stable identifiers)
  ##
  ## Description:
  ##   this function returns hgnc ids for protein coding genes symbols
  ##   however, the input file could be modified to include other genes
  ##
  ## Input:
  ##   gene_symbols: a vector of gene symbols to be checked
  ##   gene.file: a data frame 'gene_with_protein_product.txt' (see README)
  ##              file from the hgnc website
  ##
  ## Output:
  ##   A dataframe with 3 columns:
  ##     c1) 'hgnc_id': corresponding hgnc id ('-' , if no hgnc id was found)
  ##     c2) 'gene_symbol': gene symbol provided
  ##     c3) 'type': mapping type (approved_symbol,
  ##                               synonym_symbol,
  ##                               Notfound.ProteinCoding.Symbol,...)
  
  
  check.approved <- function(input.genes, database) {
    ## A function to check if the input matches an approved gene symbol
    
    return(database %>%
             dplyr::select(hgnc_id, symbol) %>%
             mutate_if(is.factor, as.character) %>%
             filter(symbol != "") %>%
             filter(!is.na(symbol)) %>%
             filter(symbol %in%
                      input.genes) %>%
             dplyr::rename(gene_symbol = symbol, hgnc_id = hgnc_id) %>%
             mutate(type = "approved_symbol"))
  }
  
  
  check.synonyms <- function(input.genes, database) {
    ## A function to check if the input symbol corresponds to an
    ## alias / synonym gene symbol
    
    return(database %>%
             dplyr::select(hgnc_id, alias_symbol) %>%
             mutate_if(is.factor,
                       as.character) %>%
             filter(alias_symbol != "") %>%
             filter(!is.na(alias_symbol)) %>%
             separate_rows(alias_symbol, sep = "\\|") %>%
             dplyr::rename(hgnc_id = hgnc_id) %>%
             mutate(gene_symbol = trimws(alias_symbol)) %>%
             dplyr::select(hgnc_id, gene_symbol) %>%
             filter(gene_symbol %in% input.genes) %>%
             mutate(type = "synonym_symbol"))
  }
  
  
  check.previous <- function(input.genes, database) {
    ## A function to check if the input symbol corresponds to a previous
    ## official gene symbol
    
    return(database %>%
             dplyr::select(hgnc_id, prev_symbol) %>%
             mutate_if(is.factor,
                       as.character) %>%
             filter(prev_symbol != "") %>%
             filter(!is.na(prev_symbol)) %>%
             separate_rows(prev_symbol, sep = "\\|") %>%
             dplyr::rename(hgnc_id = hgnc_id) %>%
             mutate(gene_symbol = trimws(prev_symbol)) %>%
             dplyr::select(hgnc_id,
                           gene_symbol) %>%
             filter(gene_symbol %in% input.genes) %>%
             mutate(type = "previous_symbol"))
  }
  
  
  check.duplicates.symbol <- function(file.to.check.symbols,
                                      duplicates.symbol) {
    ## A function to check if we have any duplicated symbols
    
    if (!length(duplicates.symbol)) {
      return(file.to.check.symbols)
      
    } else {
      
      final.nodup.symbol <- file.to.check.symbols %>%
        filter(!gene_symbol %in%
                 duplicates.symbol)
      
      duplicate.symbols.df <- data.frame(hgnc_id = rep("-",
                                                       length(duplicates.symbol)),
                                         gene_symbol = duplicates.symbol,
                                         type = "ambiguous_symbol")
      
      final.dups.symbol <- rbind(final.nodup.symbol, duplicate.symbols.df)
      
      return(final.dups.symbol)
      
    }
  }
  
  
  check.duplicates.id <- function(file.to.check.ids, duplicates.id) {
    ## A function to check if we have any duplicated hgnc ids in the
    ## resulting file
    
    if (!length(duplicates.id)) {
      return(file.to.check.ids)
    } else {
      
      final.nodup.id <- file.to.check.ids %>%
        filter(!hgnc_id %in% duplicates.id)
      
      duplicate.ids.df <- file.to.check.ids %>%
        filter(hgnc_id %in% duplicates.id) %>%
        mutate(hgnc_id = "-", type = "ambiguous_symbol")
      
      final.dups.id <- rbind(final.nodup.id, duplicate.ids.df)
      
      return(final.dups.id)
      
    }
  }
  
  # we make sure we remove any leading or trailing white space in our vector
  # with gene names
  
  genes <- trimws(gene_symbols)
  
  hgnc <- gene.file
  
  # we first check with gene symbols match the official gene symbol
  
  hgnc.approved <- check.approved(genes, hgnc)
  
  
  # we next check with gene symbols not matching the official gene symbols are
  # synonyms or alias of an approved symbol
  
  hgnc.synonyms <- check.synonyms(genes[!genes %in% hgnc.approved$gene_symbol],
                                  hgnc)
  
  
  # we next check with gene symbols not matching the official gene symbols are
  # synonyms or alias of an approved symbol
  
  hgnc.previous <- check.previous(genes[!genes %in% c(hgnc.approved$gene_symbol,
                                                      hgnc.synonyms$gene_symbol)],
                                  hgnc)
  
  
  # we next identified with genes have not been identified in the hgnc file
  # and create a data frame with 3 columns
  
  genes.not.found <- genes[!genes %in% c(hgnc.approved$gene_symbol,
                                         hgnc.synonyms$gene_symbol,
                                         hgnc.previous$gene_symbol)]
  
  
  hgnc.notfound <- data.frame(hgnc_id = rep("-", length(genes.not.found)),
                              gene_symbol = genes.not.found) %>%
    mutate(type = "Notfound.ProteinCoding.Symbol") %>%
    mutate_if(is.factor, as.character)
  
  
  # we bind all the previous data frames
  
  hgnc.all <- hgnc.approved %>%
    bind_rows(hgnc.synonyms) %>%
    bind_rows(hgnc.previous) %>%
    bind_rows(hgnc.notfound)
  
  
  # we look for any potential duplicated gene symbols or hgnc ids in our
  # dataset
  
  duplicates.symbol <- hgnc.all %>%
    filter(duplicated(gene_symbol)) %>%
    pull(gene_symbol)
  
  
  results.noduplicated.symbol <- check.duplicates.symbol(hgnc.all,
                                                         duplicates.symbol)
  
  
  duplicates.id <- results.noduplicated.symbol %>%
    filter(duplicated(hgnc_id)) %>%
    filter(hgnc_id != "-") %>%
    pull(hgnc_id)
  
  results.noduplicated.id <- check.duplicates.id(results.noduplicated.symbol,
                                                 duplicates.id)
  
  
  # final dataframe
  
  results.final <- results.noduplicated.id
  
  return(results.final)
}