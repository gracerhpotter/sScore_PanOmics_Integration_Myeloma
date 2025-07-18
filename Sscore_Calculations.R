#' @title S-Score Analysis
#' 
#' @description Application of S-score integration.
#' 
#' Citation for original publication: Nat Commun. 2013:4:2617. doi: 10.1038/ncomms3617

pacman::p_load(tidyverse,
               org.Hs.eg.db,
               AnnotationDbi,
               stringr,
               dplyr,
               readr)

prefix__id_proteome <- "pr."
prefix__id_rna <- "rna."
prefix__id_phosphoproteome <- "ph."
prefix__id_metabolome <- "met."

pcsf_pval_cutoff <- 0.05
my_pval <- pcsf_pval_cutoff
label_pval <- pcsf_pval_cutoff

# Function for computing S-scores for genes
compute_sscore_for_genes <- function(x) {
  x %>%
    dplyr::mutate(
      uniprot_id = stringr::str_squish(uniprot_id),
      comb_weighted_zi = rowSums(dplyr::select(., dplyr::starts_with("weighted_zi_")), na.rm = TRUE),
      comb_wk = sqrt(rowSums(dplyr::select(., dplyr::starts_with("wk_"))^2, na.rm = TRUE))
    ) %>%
    dplyr::mutate(
      sscore = comb_weighted_zi / comb_wk,
      sscore_pval = stats::pnorm(abs(sscore), lower.tail = FALSE) * 2,
      sscore_adj_pval = stats::p.adjust(sscore_pval, method = 'BH')
    ) %>%
    dplyr::mutate(across(where(is.character), ~ na_if(., "NA"))) %>%
    as.data.frame() %>%
    dplyr::select(
      uniprot_id, 
      dplyr::starts_with("feature_id"), 
      dplyr::starts_with("gene_symbol"),
      sscore, 
      sscore_pval, 
      sscore_adj_pval,
      dplyr::starts_with("comb_"),
      dplyr::starts_with("logfc"), 
      dplyr::starts_with("weighted"),
      dplyr::starts_with("wk")
    )
}

# Function for computing S-scores for metabolites
compute_sscore_for_metabols <- function(x) {
  x %>%
    dplyr::mutate(
      chebi_id = stringr::str_squish(chebi_id),
      comb_weighted_zi = rowSums(dplyr::select(., dplyr::starts_with("weighted_zi_")), na.rm = TRUE),
      comb_wk = sqrt(rowSums(dplyr::select(., dplyr::starts_with("wk_"))^2, na.rm = TRUE)),
      sscore = comb_weighted_zi / comb_wk,
      sscore_pval = stats::pnorm(abs(sscore), lower.tail = FALSE) * 2,
      sscore_adj_pval = stats::p.adjust(sscore_pval, method = 'BH'),
      feature_id_metabolome = chebi_id,
      uniprot_id = "NA",
      gene_symbol = "NA",
    ) %>%
    dplyr::mutate(across(where(is.character), ~ na_if(., "NA"))) %>%
    as.data.frame() %>%
    dplyr::select(
      uniprot_id,
      chebi_id,
      dplyr::starts_with("feature_id"),
      dplyr::starts_with("gene_symbol"),
      sscore,
      sscore_pval,
      sscore_adj_pval,
      dplyr::starts_with("comb_"),
      dplyr::starts_with("logfc"),
      dplyr::starts_with("weighted"),
      dplyr::starts_with("wk")
    )
}

# Define the function to compute weighted Z-scores and S-scores
compute_weighted_zi_list <- function(data_list = data_list) {
  
  message(paste("There are", length(data_list), "elements in data_list") )
  message(" ")
  message("*** Input data format ***")
  message(" -- The input 'data_list' is a list of 'dataframes' with 3 columns: uniprot_id, feature_id and logfc")
  message(" -- The code will combine all elements of 'data_list' using 'uniprot_id'.")
  message(" -- If present, S-score for metabolite dataset will be computed separately and then combined with genes")
  message(" -- Use 'my_combine_genes_and_metabols_using_sscore' function if metabolomics data is present.")
  message(" -- Use 'my_combine_genes_using_sscore' function if metabolomics data is NOT present.")
  message(" -- Example 'feature_id' for proteomics: Q8BH50_ARK2N")
  message(" -- Example 'feature_id' for phosphoproteomics & other PTMs: Q8BH50_ARK2N_T74")
  message(" ")
  message("*** Citation ***")
  message("Citation for original publication: Nat Commun. 2013:4:2617. doi: 10.1038/ncomms3617" )
  
  # Define the biological importance weights
  bio_importance_weights <- list(protein = 0.5, 
                                 phosphosite = 0.5, 
                                 rna = 0.3)
  
  # Calculate the number of features for each dataset
  n_protein <- nrow(data_list$proteome)
  n_phosphosite <- nrow(data_list$phosphoproteome)
  n_rna <- nrow(data_list$rna)
  
  # Calculate weights based on the number of features and biological importance
  protein_weight <- bio_importance_weights$protein / sqrt(n_protein)
  phosphosite_weight <- bio_importance_weights$phosphosite / sqrt(n_phosphosite)
  rna_weight <- bio_importance_weights$rna / sqrt(n_rna)
  
  # Internal function to compute weighted Z-scores for each dataset
  compute_weighted_zi <- function(x, dataset_type) {
    wk <- switch(dataset_type,
                 "proteome" = protein_weight,
                 "phosphoproteome" = phosphosite_weight,
                 "rna" = rna_weight)
    
    x %>%
      dplyr::mutate(wk = wk,
                    zi = scale(logfc, center = TRUE, scale = TRUE)[,1],
                    weighted_zi = wk * zi) %>%
      dplyr::select(logfc, wk, weighted_zi, everything()) %>%
      dplyr::mutate(across(where(is.character), ~ na_if(., "NA"))) %>%
      as.data.frame()
  }
  
  for (i in 1:length(data_list)) {
    
    # Extract data frame from the list
    data_s <- data_list[[i]]
    
    # Extract element name from the list
    element_name <- names(data_list)[i]
    
    # Compute weighted Z-scores using the new approach
    data_s1 <- compute_weighted_zi(data_s, element_name) %>%
      dplyr::rename_with(~ ifelse(.x %in% c("uniprot_id", "gene_symbol"), .x, paste0(.x, "_", element_name)))
    
    # Save the computed data frame back into the list
    data_list[[i]] <- data_s1
    
  }
  
  # Assign the modified list to the global environment
  assign("data_list_mod", data_list, envir = .GlobalEnv)
  
  return(data_list)
}

# final combine using S-score
my_combine_genes_and_metabols_using_sscore <- function(data_list = data_list, 
                                                       pcsf_pval_cutoff = pcsf_pval_cutoff, 
                                                       output_name = output_name) {
  
  compute_weighted_zi_list(data_list)
  
  # Combine all data frames in the GENES list using compute_sscore function
  elements_to_use <- names(data_list_mod)[names(data_list_mod) != "metabolome"]
  sscore_genes <- data_list_mod[elements_to_use] %>% 
    purrr::reduce(full_join, by = "gene_symbol") %>% # combining with SYMBOL will expand the rows because multiple UNIPROT per gene
    compute_sscore_for_genes() %>%
    dplyr::arrange(sscore_adj_pval)
  
  assign("sscore_genes", sscore_genes, envir = .GlobalEnv)
  
  # Bind rows
  sscore_combined_genes_metabols <- sscore_genes %>%
    dplyr::mutate(across(where(is.character), ~ na_if(., "NA"))) %>%
    as.data.frame() %>%
    dplyr::select(uniprot_id,
                  dplyr::starts_with("feature_id"), 
                  dplyr::starts_with("gene_symbol"), 
                  sscore, 
                  sscore_pval, 
                  sscore_adj_pval, 
                  dplyr::starts_with("comb_"), 
                  dplyr::starts_with("logfc"), 
                  dplyr::starts_with("weighted"), 
                  dplyr::starts_with("wk") ) %>%
    dplyr::mutate(output_label = output_name) 
  
  columns_to_process <- sscore_combined_genes_metabols %>% dplyr::select(starts_with("feature_id_")) %>% names()
  
  sscore_combined_genes_metabols <- sscore_combined_genes_metabols %>%
    as.data.frame() %>%
    dplyr::mutate(across(where(is.character), ~ na_if(., "NA"))) %>%
    dplyr::mutate(across(where(is.character), ~ na_if(., "NULL"))) %>%
    as.data.frame() %>%
    dplyr::mutate(dplyr::across(dplyr::all_of(columns_to_process), ~stringr::str_replace(., ".*_", ""))) %>%
    dplyr::mutate(dplyr::across(dplyr::contains("feature_id_"), 
                                ~{
                                  column_name <- as.character(cur_column())
                                  matching_string <- str_extract(column_name, "_(.*)")
                                  
                                  if (!is.na(matching_string)) {
                                    prefix_variable <- paste0("prefix_", tolower(gsub("\\d", "", matching_string)))
                                    paste(get(prefix_variable), ., sep = "")
                                  } else {
                                    .  # For columns without matching conditions, keep the original value
                                  }
                                }
    )) %>%
    dplyr::mutate(across(where(is.character), ~ na_if(., "NA"))) %>%
    dplyr::mutate(across(where(is.character), ~ na_if(., "NULL"))) %>%
    as.data.frame() %>%
    dplyr::rowwise() %>%
    dplyr::mutate(sscore_label = paste(c(gene_symbol, dplyr::across(dplyr::all_of(columns_to_process), as.character)), collapse = "_"),
                  feature_id = uniprot_id) %>%
    dplyr::ungroup() %>%
    as.data.frame() %>%
    dplyr::select(sscore_label, feature_id, gene_symbol, sscore, sscore_pval, sscore_adj_pval, everything())
  
  assign("sscore_combined_genes_metabols", sscore_combined_genes_metabols, envir = .GlobalEnv)
  
  # write file
  folder_path <- paste0(getwd(), "/sscore_output/")
  
  # Create the folder if it doesn't exist
  if (!file.exists(folder_path)) {
    dir.create(folder_path, recursive = TRUE)
  }
  
  readr::write_delim(sscore_combined_genes_metabols, paste0(folder_path, "sscore_full_", output_name, ".txt"), delim = "\t")
  
  # OUTPUT FORMATTED FOR PCSF
  out_for_pcsf <- sscore_combined_genes_metabols %>%
    dplyr::filter(sscore_adj_pval <= pcsf_pval_cutoff) %>%
    dplyr::mutate(logfc = sscore * -log10(sscore_adj_pval)) %>% # Remove S-score ties
    dplyr::arrange(desc(abs(logfc))) %>%
    dplyr::group_by(sscore_label) %>%
    dplyr::slice_max(order_by = abs(logfc), n = 1) %>%
    dplyr::ungroup() %>%
    dplyr::select(gene_symbol, sscore_label, logfc, sscore, sscore_adj_pval) %>%
    dplyr::distinct(gene_symbol, .keep_all = TRUE) %>%
    readr::write_delim(paste0(folder_path, "sscore_pcsf_prize_pval", pcsf_pval_cutoff, "_",output_name, ".txt"), delim = "\t")
}

################################################################################

setwd("")

# LOAD DATA
dataList_H9 <- list(proteome = readr::read_delim("sscoreInputs/Prot_H9.txt", delim = "\t"),
                    phosphoproteome = readr::read_delim("sscoreInputs/Phospho_H9.txt", delim = "\t"),
                    rna = readr::read_delim("sscoreInputs/RNA_biomaRt_H9.txt", delim = "\t"))
dataList_JI <- list(proteome = readr::read_delim("sscoreInputs/Prot_JI.txt", delim = "\t"),
                    phosphoproteome = readr::read_delim("sscoreInputs/Phospho_JI.txt", delim = "\t"),
                    rna = readr::read_delim("sscoreInputs/RNA_biomaRt_JI.txt", delim = "\t"))
dataList_KM <- list(proteome = readr::read_delim("sscoreInputs/Prot_KM.txt", delim = "\t"),
                    phosphoproteome = readr::read_delim("sscoreInputs/Phospho_KM.txt", delim = "\t"),
                    rna = readr::read_delim("sscoreInputs/RNA_biomaRt_KM.txt", delim = "\t"))

# RUN
my_combine_genes_and_metabols_using_sscore(data_list = dataList_H9, 
                                           pcsf_pval_cutoff = 0.05, 
                                           output_name = "H9_diffWeights") 

my_combine_genes_and_metabols_using_sscore(data_list = dataList_JI, 
                                           pcsf_pval_cutoff = 0.05, 
                                           output_name = "JI_diffWeights") 

my_combine_genes_and_metabols_using_sscore(data_list = dataList_KM, 
                                           pcsf_pval_cutoff = 0.05, 
                                           output_name = "KM_diffWeights") 

