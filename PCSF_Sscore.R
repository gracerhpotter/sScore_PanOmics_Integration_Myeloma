
# S-SCORE ######################################################################
## FUNCTIONS -------------------------------------------------------------------
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

# Now we can define the function to compute weighted Z-scores and S-scores
compute_weighted_zi_list <- function(data_list = data_list) {
  
  message(paste("There are", length(data_list), "elements in data_list") )
  message(" ")
  message("*** Input data format ***")
  message(" -- The input 'data_list' is a list of 'dataframes' with 3 columns: uniprot_id/gene_symbol, feature_id and logfc")
  message(" -- The code will combine all elements of 'data_list' using 'uniprot_id/gene_symbol'.")
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
  sscore_combined_genes <- sscore_genes %>%
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
  
  columns_to_process <- sscore_combined_genes %>% dplyr::select(starts_with("feature_id_")) %>% names()
  
  sscore_combined_genes <- sscore_combined_genes %>%
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
  
  # write file
  folder_path <- paste0(getwd(), "/sscore_output/")
  
  # Create the folder if it doesn't exist
  if (!file.exists(folder_path)) {
    dir.create(folder_path, recursive = TRUE)
  }
  
  readr::write_delim(sscore_combined_genes, paste0(folder_path, "sscore_full_", output_name, ".txt"), delim = "\t")
  
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

## RUN -------------------------------------------------------------------------
# LOAD PACKAGES
library("pacman")
pacman::p_load(tidyverse,
               org.Hs.eg.db,
               AnnotationDbi,
               VennDiagram,
               readr,
               dplyr,
               stringr)

prefix__id_proteome <- "pr."
prefix__id_rna <- "rna."
prefix__id_phosphoproteome <- "ph."

pcsf_pval_cutoff <- 0.05
my_pval <- pcsf_pval_cutoff
label_pval <- pcsf_pval_cutoff
my_sscore <- 3.7

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

# PCSF #########################################################################
## FUNCTIONS -------------------------------------------------------------------

# Prize-Collecting Steiner Forest (PCSF) Graph Optimization Approach
# PLoS Comput Biol. 2017 Jul; 13(7): e1005694. DOI: 10.1371/journal.pcbi.1005694.

# PCSF NETWORK ANALYSIS ########################################################

## MAKE PCSF NETWORK -----------------------------------------------------------
#' @title PCSF Network
#' 
#' @description
#' Create PCSF network using PCSF_rand function from PCSF package for use in createing PCSF
#' network visualization.
#' 
#' @param pcsf_input_data dataframe input - must contain gene_symbol column and value column (logFC or S-score)
#' @param pcsf_nval numeric input -  number of runs with random noise added edge costs
#' 
#' @return PCSF network, PCSF igraph list object
#'

makePCSFNetwork <- function(pcsf_input_data = pcsf_input_data,
                            pcsf_nval = 10,
                            mu = 0.005){
  
  pcsf_ppi <- readr::read_rds("pcsf_ppi.rds")
  
  names(pcsf_input_data)[4] <- 'value'
  
  # input prizes
  pcsf_input_prizes = pcsf_input_data %>% 
    dplyr::select(gene_symbol, value) %>% 
    dplyr::mutate(gene_symbol = stringr::str_squish(gene_symbol))
  
  pcsf_input_prizes$value = abs(pcsf_input_prizes$value) # PCSF takes only positive values
  pcsf_input_prizes = tibble::deframe(pcsf_input_prizes) 
  pcsf_input_prizes_df = tibble::enframe(pcsf_input_prizes)
  
  set.seed(123)
  pcsf_net = PCSF::PCSF_rand(ppi = pcsf_ppi, 
                             terminal = pcsf_input_prizes, 
                             n = pcsf_nval,
                             mu = mu) # higher mu = smaller the number of Steiners
  
  return(pcsf_net)
}

## INTERACTION NETWORK ---------------------------------------------------------
#' @title PCSF Interaction Network
#'
#' @description
#' Generates an interactive PCSF network where nodes are genes and edges are interactions.
#' 
#' @param pcsf_net PCSF network output from makePCSFNetwork()
#' @param pcsf_input_data dataframe input
#'
#' @return visNetwork::visIgraph object

### PERFORM PCSF NETWORK ANALYSIS
PCSFVisNodes = function(pcsf_net = pcsf_net,
                        pcsf_input_data = pcsf_input,
                        name = "") {
  
  my_interactome <- readr::read_rds("pcsf_interactome.rds")
  names(pcsf_input_data)[4] <- 'value'
  
  # EDGES
  my_pcsf_net_edges = as.data.frame(igraph::get.edgelist(pcsf_net)) %>%
    dplyr::mutate(COMBINATION = paste0(pmin(V1, V2), pmax(V1, V2))) %>%
    dplyr::left_join(., my_interactome, by = "COMBINATION") %>%
    dplyr::mutate(weights = WEIGHT) %>%
    dplyr::mutate(weights = tidyr::replace_na(weights, min(weights, na.rm = TRUE))) %>%
    dplyr::mutate(label = INTERACTION_TYPE) %>%
    dplyr::mutate(label = replace(label, label == "interacts-with", "")) %>%
    dplyr::select(V1, V2, weights, label)
  
  # NODES
  my_pcsf_net_nodes = data.frame(gene_symbol = c(my_pcsf_net_edges$V1, my_pcsf_net_edges$V2))
  my_pcsf_net_nodes = my_pcsf_net_nodes %>% 
    dplyr::distinct(gene_symbol) %>%
    dplyr::left_join(., pcsf_input_data, by = "gene_symbol") %>%
    dplyr::distinct(gene_symbol, .keep_all = TRUE) %>%
    dplyr::mutate(is_steiner = case_when(value == "NA" ~ "Yes", 
                                         value != "NA" ~ "No")) %>% 
    dplyr::mutate(is_steiner = tidyr::replace_na(is_steiner, "Yes")) %>% # if no logfc then its a Steiner
    dplyr::mutate(regulation = case_when(value >= 0 ~ "Up",
                                         value <= -0 ~ "Down")) %>%
    dplyr::mutate(regulation = tidyr::replace_na(regulation, "None")) %>%
    dplyr::mutate(value = tidyr::replace_na(value, 0)) %>%
    dplyr::mutate(abs_value = log2(abs(value))) 
  
  steiner_value = min(my_pcsf_net_nodes$abs_value[is.finite(my_pcsf_net_nodes$abs_value)], na.rm = TRUE) 
  my_pcsf_net_nodes = my_pcsf_net_nodes %>% 
    dplyr::mutate_all(~ifelse(is.na(.x) | is.nan(.x) | .x == -Inf, steiner_value, .x))
  
  file <- paste0("~/Downloads/", name, "_pcsf_nodeEdgeTable.xlsx")
  wb <- xlsx::createWorkbook()
  sheet_node <- xlsx::createSheet(wb, sheetName = "Nodes")
  sheet_edge <- xlsx::createSheet(wb, sheetName = "Edges")
  xlsx::addDataFrame(my_pcsf_net_nodes, sheet = sheet_node)
  xlsx::addDataFrame(my_pcsf_net_edges, sheet = sheet_edge)
  xlsx::saveWorkbook(wb, file = file)
  
  # Define igraph nodes and edge aesthetics. Use igraph-compatible names
  my_pcsf_net_nodes = my_pcsf_net_nodes %>%
    dplyr::mutate(color = case_when(regulation == "Up" ~ "#F14E2B",
                                    regulation == "Down" ~ "#3F5DF4",
                                    regulation == "None" ~ "#E7E7E8")) %>%
    dplyr::mutate(label = gene_symbol)
    # dplyr::mutate(label.cex = scales::rescale(abs_value, c(0.3, 0.9)))
  
  # create igraph object
  my_pcsf_net_vis = igraph::graph_from_data_frame(d = my_pcsf_net_edges, 
                                                  vertices = my_pcsf_net_nodes, 
                                                  directed = FALSE)
  igraph::write_graph(graph = my_pcsf_net_vis,
                      file = paste0("~/Downloads/", name, "_pcsf_nodenet.graphml"),
                      format = "graphml")
  # plot
  set.seed(123)
  visIgraph_obj = visNetwork::visIgraph(igraph = my_pcsf_net_vis,
                                        layout = "layout_nicely",
                                        physics = TRUE,
                                        smooth = TRUE, 
                                        idToLabel = FALSE) %>%
    visNetwork::visNodes(shape = "circle")  %>%
    visNetwork::visEdges(smooth = TRUE,
                         value = 2) %>%
    visNetwork::visOptions(height = "1000px",
                           width = "1000px",
                           highlightNearest = TRUE,
                           nodesIdSelection = TRUE ) %>%
    visNetwork::visInteraction(navigationButtons = TRUE, 
                               dragNodes = TRUE,
                               dragView = TRUE, 
                               zoomView = TRUE) 
  
  visNetwork::visSave(visIgraph_obj, file = paste0("~/Downloads/", name, "_pcsf_nodenet.html"), selfcontained = TRUE, background = "white")
}
  
## INFLUENTIAL NETWORK ---------------------------------------------------------
#' @title PCSF Influential Node Network
#' 
#' @description
#' Generates an interactive PCSF network where nodes are genes and they are colored
#' and sized based on measure of "influence" determined by the `influential` package.
#' 
#' @param pcsf_net PCSF network object
#' 
#' @return visNetwork::visIgraph object
#'

PCSFVisInfluential = function(pcsf_net = pcsf_net,
                              name = "") {
  ### PERFORM INFLUENTIAL ANALYSIS ON PCSF NETWORK
  my_graph = pcsf_net
  my_graph_vertices = igraph::V(pcsf_net) 
  
  # Compute influential nodes
  set.seed(123)
  my_graph_ivi = influential::ivi(graph = my_graph, 
                                  vertices = my_graph_vertices, 
                                  weights = NULL, 
                                  directed = FALSE, 
                                  mode = "all",
                                  loops = TRUE, 
                                  d = 3, 
                                  scale = "range" )
  
  my_graph_ivi_df = data.frame(my_graph_ivi) %>% 
    dplyr::mutate(genes = rownames(.)) %>%
    dplyr::rename(influence_score = my_graph_ivi) %>%
    dplyr::select(genes, influence_score) %>%
    dplyr::arrange(desc(influence_score))
  
  readr::write_tsv(my_graph_ivi_df, paste0("~/Downloads/", name, "_pcsf_influential.tsv"))
  
  # assign("PCSF_influential", my_graph_ivi_df, envir = .GlobalEnv)
  
  my_graph_ivi_edges = as.data.frame(igraph::get.edgelist(my_graph))
  my_graph_ivi_nodes = my_graph_ivi_df
  
  # create igraph object
  my_igraph_ivi_vis = igraph::graph_from_data_frame(d = my_graph_ivi_edges, 
                                                    vertices = my_graph_ivi_nodes, 
                                                    directed = FALSE) 
  
  # pass node and edge information
  igraph::V(my_igraph_ivi_vis)$label <- igraph::as_ids(igraph::V(my_igraph_ivi_vis))
  igraph::V(my_igraph_ivi_vis)$color <- colourvalues::colour_values(igraph::V(my_igraph_ivi_vis)$influence_score, palette = "spectral", include_alpha = FALSE)
  
  igraph::V(my_igraph_ivi_vis)$label.cex <- igraph::V(my_igraph_ivi_vis)$influence_score
  igraph::V(my_igraph_ivi_vis)$label.cex <- scales::rescale(igraph::V(my_igraph_ivi_vis)$label.cex, c(0.5, 1.5))
  igraph::V(my_igraph_ivi_vis)$size <- igraph::V(my_igraph_ivi_vis)$influence_score
  igraph::V(my_igraph_ivi_vis)$size <- scales::rescale(igraph::V(my_igraph_ivi_vis)$size, c(8, 24))
  
  igraph::write_graph(graph = my_igraph_ivi_vis,
                      file = paste0("~/Downloads/", name, "_pcsf_influential.graphml"),
                      format = "graphml")
  
  # plot
  visIgraph_obj = visNetwork::visIgraph(igraph = my_igraph_ivi_vis,
                                        layout = "layout_nicely",
                                        physics = TRUE,
                                        smooth = TRUE, 
                                        idToLabel = FALSE) %>%
    visNetwork::visNodes() %>%
    visNetwork::visEdges() %>%
    visNetwork::visOptions(height = "1000px",
                           width = "1000px",
                           highlightNearest = TRUE,
                           nodesIdSelection = TRUE ) %>%
    visNetwork::visInteraction(navigationButtons = TRUE, 
                               dragNodes = TRUE,
                               dragView = TRUE, 
                               zoomView = TRUE) 
  visNetwork::visSave(visIgraph_obj, file = paste0("~/Downloads/", name, "_pcsf_influential.html"), selfcontained = TRUE, background = "white")
  
}

# PCSF ENRICHMENT ##############################################################

## CLUSTER PROFILER ------------------------------------------------------------
#' @title PCSF Run ClusterProfiler
#'
#' @description
#' Run clusterProfiler and format output.
#' 
#' @param clusters igraph::cluster_louvain(subnet)
#' @param subnet PCSF network
#' @param gmt name of GMT database file.
#' 
#' @return enrichment_result_complete

runClusterProfiler <- function(clusters,
                               subnet,
                               gmt,
                               name, 
                               wb) {
  
  # gmt = "HUMAN_KEGG"
  enrichment_result = as.list(1:length(clusters))
  enrichment_result_complete = as.list(1:length(clusters))
  enrichment_table = as.list(1:length(clusters))
  
  # GENERATE GMT REFERENCES
  my_geneset <- clusterProfiler::read.gmt("Human_Reactome_January_04_2025_symbol.gmt")
  colnames(my_geneset) <- c("pathway", "feature_ids")
  
  my_geneset = dplyr::mutate(my_geneset, term = pathway) %>% dplyr::mutate(feature = stringr::str_squish(feature_ids))
  term2features = my_geneset %>% dplyr::select(term, feature) #TERM2GENE
  term2pathway = my_geneset %>% dplyr::select(term, pathway) #TERM2NAME
  length(unique(my_geneset$feature))
  
  set.seed(123)
  enriched_dataList <- list()
  for (a in 1:length(clusters)) {
    
    # INDICATE CLUSTER NUMBER
    cluster_number <- a
    print(paste0("Cluster Number: ", a))
    
    assign("clusters", clusters, envir = .GlobalEnv)
    
    enrichment_input <- clusters[[a]]
    background_genes <- unique(igraph::V(subnet)$name) # use network as background
    
    tryCatch({
      enriched_pathways <- clusterProfiler::enricher(gene = enrichment_input,
                                                     universe = background_genes,
                                                     minGSSize = 2,
                                                     pAdjustMethod = "fdr",
                                                     pvalueCutoff = 0.05,
                                                     TERM2GENE = term2features,
                                                     TERM2NAME = term2pathway)
      enriched_dataList[[a]] <- enriched_pathways
      
      # EMPTY DATAFRAME TO STORE OUTPUT
      enriched_pathways_df_a = data.frame(cluster = "NULL", 
                                          term = "NULL", 
                                          adj_pval = 1, 
                                          intersection = "NULL")
      # DATAFRAME WITH ENRICHER RESULTS
      enriched_pathways_df_b = base::data.frame(enriched_pathways@result) %>%
        tibble::remove_rownames() %>%
        dplyr::select(term = Description, 
                      adj_pval = p.adjust, 
                      intersection = geneID) %>%
        dplyr::mutate(intersection = stringr::str_replace_all(intersection, "/", ", "))
      
      # BIND DATAFRAMES
      enriched_pathways_df = dplyr::bind_rows(enriched_pathways_df_a, enriched_pathways_df_b) %>%
        dplyr::mutate(cluster = paste0("Cluster_", cluster_number) )
      
      # DATA TABLE OF TOP RESULTS
      res_table_top15 = enriched_pathways_df %>% 
        dplyr::top_n(n = 15, wt = -adj_pval) %>%
        dplyr::filter(adj_pval <= 0.05) %>%
        dplyr::select(term, adj_pval, intersection) %>%
        dplyr::arrange(adj_pval)
      
      enrich = "<!DOCTYPE html> <html> <head> <style>\n      table {font-family: arial, sans-serif; font-size: 10px; border-collapse: collapse;width: 100%;} td,\n      th { border: 1px solid #dddddd; text-align: center; padding: 5px;}\n      tr:nth-child(even) {background-color: #dddddd;}\n      </style> </head> <body>\n      <table> <tr>  <th>Term</th> <th>Adjusted P-value</th> <th>Intersection</th> </tr>"
      
      for (i in 1:nrow(res_table_top15)) {
        enrich = paste0(enrich, " <tr>")
        for (j in 1:ncol(res_table_top15)) {
          enrich = paste0(enrich, "<td>", res_table_top15[i, 
                                                          j], "</td>")
        }
        enrich = paste0(enrich, "</tr> ")
      }
      
      enrich = paste0(enrich, "</table> </body> </html>")
      enrichment_result[[a]] = enrich
      enrichment_result_complete[[a]] = enriched_pathways_df
      enrichment_table[[a]] = enriched_pathways_df_b
      
      # Further processing of enriched_pathways if needed
    }, error = function(e) {
      warning(paste("Error in enricher for cluster ", a, ": ", conditionMessage(e)))
    })
  }
  
  df <- NULL
  enrichment_summary <- data.frame()
  for (i in 1:length(enriched_dataList)){
    try(df <- enriched_dataList[[i]]@result)
    if (!is.null(df)){
      df[,1] <- rep(i, nrow(df))
      enrichment_summary <- rbind(enrichment_summary, df)
    }
    
    df <- NULL
  }
  
  sheet <- xlsx::createSheet(wb, sheetName = "Enrichment_Summary")
  xlsx::addDataFrame(enrichment_summary, sheet = sheet)
  
  return(list(enrichment_result, enrichment_result_complete, enrichment_table))
}

## ENRICHMENT ------------------------------------------------------------------
#' @title PCSF Enrichment
#'
#' @description
#' Call gprofiler and map pathways to PCSF modules.
#' 
#' @param subnet pcsf network
#' @param gmt GMT database name
#' 
#' @return output = list(subnet, enrichment_tab)

PCSFModuleEnrichments <- function(subnet,
                                  gmt,
                                  name, 
                                  wb) {
  
  # Perform clustering on the PCSF network
  set.seed(123)
  clusters = igraph::cluster_louvain(subnet)
  # clusters = igraph::cluster_fast_greedy(subnet) 
  # clusters = igraph::cluster_edge_betweenness(subnet) 
  
  cluster_df <- data.frame("Cluster_1", paste0(clusters[[1]], collapse = "; "))
  
  for (i in 2:length(clusters)){
    cluster_df <- rbind(cluster_df, c(paste0("Cluster_", i), paste0(clusters[[i]], collapse = "; ")))
  }
  
  colnames(cluster_df) <- c("Cluster", "Genes")
  
  sheet_cluster <- xlsx::createSheet(wb, sheetName = "Clusters")
  xlsx::addDataFrame(cluster_df, sheet = sheet_cluster)
  
  enrich = runClusterProfiler(clusters,
                              subnet,
                              gmt,
                              name,
                              wb)
  
  enrichment = enrich[[1]]
  enrichment_complete = enrich[[2]]
  
  novals = which(unlist(sapply(enrich[[2]], function(x) is.null(dim(x)))))
  
  if (length(novals) > 0){
    enrichment_complete = enrichment_complete[-novals]
  } 
  
  if (length(enrichment_complete) == 0){
    return(NULL)
  }
  
  enrichment_tab = do.call(rbind, lapply(c(1:length(enrichment_complete)), function(x) data.frame(Cluster = x, enrichment_complete[[x]])))
  
  igraph::V(subnet)$group = clusters$membership
  igraph::V(subnet)$title = paste0("Cluster ", clusters$membership, ": Enrichment analysis")
  
  for (i in 1:length(igraph::V(subnet))) {
    igraph::V(subnet)$title[i] = paste0(igraph::V(subnet)$title[i], enrichment[[igraph::V(subnet)$group[i]]])
  }
  
  class(subnet) = c("PCSFe", "igraph")
  output = list(subnet, enrichment_tab)
  names(output) = c("subnet", "enrichment")
  
  return(output)
}

## ENRICHED NETWORK ------------------------------------------------------------
#' @title Run PCSF Enrichment
#' 
#' @description
#' Call PCSF Enrichment modules.
#' 
#' @param pcsf_net PCSF network
#' @param gmt GMT database name
#' 
#' @return Enrichment results, list enrichment and subnet.
#' 

pcsfRunEnrichment <- function(pcsf_net,
                              gmt,
                              name = ""){
  ### PERFORM PATHWAY ENRICHMENT ON PCSF NETWORK
  ## This step performs clustering and then calls enrichment on each cluster
  file <- paste0("~/Downloads/", name, "_pcsf_enrichedPathways.xlsx")
  wb <- xlsx::createWorkbook()
  
  subnet = pcsf_net
  set.seed(123)
  pcsf_enrich_pathway = PCSFModuleEnrichments(subnet,
                                              gmt,
                                              name, 
                                              wb)
  
  intersections <- pcsfEnrichedTable(pcsf_enrich_pathway, name)
  
  sheet_intersect <- xlsx::createSheet(wb, sheetName = "Enrichment_Intersections")
  xlsx::addDataFrame(intersections, sheet = sheet_intersect)
  xlsx::saveWorkbook(wb, file = file)
  
  return(pcsf_enrich_pathway)
}

#-------------------------------------------------------------------------------
#' @title PCSF Enrichment Contracted Network
#' 
#' @description
#' Visualize contracted network of PCSF enrichment results.
#' 
#' @param pcsf_enrich_pathway enrichment results output from pcsfRunEnrichment()
#' 
#' @return visNetwork::visIgraph object
#' 

pcsfEnrichedContracted <- function(pcsf_enrich_pathway){
  
  # Contract the modules
  my_g = pcsf_enrich_pathway$subnet
  contracted_graph = igraph::contract(my_g, igraph::V(my_g)$group)
  contracted_graph = igraph::simplify(contracted_graph)
  
  number_of_clusters = seq(1:length(contracted_graph))
  new_module_names = paste("Module-", number_of_clusters, sep = "")
  igraph::V(contracted_graph)$name = new_module_names
  igraph::V(contracted_graph)$group = new_module_names
  
  unique_clusters = unique(igraph::V(contracted_graph)$name)
  cluster_gene_counts = sapply(unique_clusters, length)
  igraph::V(contracted_graph)$size = cluster_gene_counts
  igraph::V(contracted_graph)$color = colourvalues::colour_values(igraph::V(contracted_graph)$name, palette = "spectral", include_alpha = FALSE)

  igraph::write_graph(graph = contracted_graph,
                      file = paste0("~/Downloads/pcsf_enriched_contracted_", Sys.Date(), ".graphml"),
                      format = "graphml")
  
  visIgraph_obj = visNetwork::visIgraph(igraph = contracted_graph,
                                        layout = "layout_nicely", #layout_nicely
                                        physics = TRUE,
                                        smooth = TRUE, 
                                        idToLabel = FALSE) %>%
    visNetwork::visNodes(font = list(size = 15),
                         shape = "circle") %>%
    visNetwork::visEdges(smooth = TRUE,
                         value = 3) %>%
    visNetwork::visOptions(height = "1000px",
                           width = "1000px",
                           highlightNearest = TRUE,
                           nodesIdSelection = TRUE,
                           selectedBy = "group") %>%
    visNetwork::visInteraction(navigationButtons = TRUE, 
                               dragNodes = TRUE,
                               dragView = TRUE, 
                               zoomView = TRUE) 
  
  return(visIgraph_obj)
}

#-------------------------------------------------------------------------------
#' @title PCSF Enriched Subnet Network
#' 
#' @description
#' Visualize subnet of PCSF enrichment results.
#'
#' @param pcsf_enrich_pathway enrichment results output from pcsfRunEnrichment()
#' 
#' @return visNetwork::visIgraph object
#'

pcsfEnrichedSubnet <- function(pcsf_enrich_pathway,
                               name = ""){
  # Create visNetwork object for nice visualization
  visIgraph_obj = pcsf_enrich_pathway$subnet
  
  igraph::V(visIgraph_obj)$cluster = gsub(":.*", "", igraph::V(visIgraph_obj)$title)
  igraph::V(visIgraph_obj)$color <- colourvalues::colour_values(igraph::V(visIgraph_obj)$cluster, palette = "spectral", include_alpha = FALSE)
  igraph::write_graph(graph = visIgraph_obj,
                      file = paste0("~/Downloads/", name, "_pcsf_enrichedSubnet.graphml"),
                      format = "graphml")
  
  visIgraph_obj_print = visNetwork::visIgraph(igraph = visIgraph_obj,
                                        layout = "layout_nicely", #layout_nicely
                                        physics = TRUE,
                                        smooth = TRUE, 
                                        idToLabel = FALSE) %>%
    visNetwork::visNodes(font = list(size = 30)) %>%
    visNetwork::visEdges(smooth = TRUE,
                         value = 3) %>%
    visNetwork::visOptions(height = "1000px",
                           width = "1000px",
                           highlightNearest = TRUE,
                           nodesIdSelection = TRUE,
                           selectedBy = "group") %>%
    visNetwork::visInteraction(navigationButtons = TRUE, 
                               dragNodes = TRUE,
                               dragView = TRUE, 
                               zoomView = TRUE)
  
  visNetwork::visSave(visIgraph_obj_print, file = paste0("~/Downloads/", name, "_pcsf_enrichedSubnet_LABEL.html"), selfcontained = TRUE, background = "white")
  
  visIgraph_obj_print = visNetwork::visIgraph(igraph = visIgraph_obj,
                                        layout = "layout_nicely", #layout_nicely
                                        physics = TRUE,
                                        smooth = TRUE, 
                                        idToLabel = FALSE) %>%
    visNetwork::visNodes(font = list(size = 0)) %>%
    visNetwork::visEdges(smooth = TRUE,
                         value = 3) %>%
    visNetwork::visOptions(height = "1000px",
                           width = "1000px",
                           highlightNearest = TRUE,
                           nodesIdSelection = TRUE,
                           selectedBy = "group") %>%
    visNetwork::visInteraction(navigationButtons = TRUE, 
                               dragNodes = TRUE,
                               dragView = TRUE, 
                               zoomView = TRUE)
  
  visNetwork::visSave(visIgraph_obj_print, file = paste0("~/Downloads/", name, "_pcsf_enrichedSubnet.html"), selfcontained = TRUE, background = "white")
}

## RUN -------------------------------------------------------------------------
# WITHOUT INTERACTOME
H9 <- read.delim("PCSF/sscore_pcsf_prize_pval0.05_H9_diffWeights.txt", sep = "\t")
JI <- read.delim("PCSF/sscore_pcsf_prize_pval0.05_JI_diffWeights.txt", sep = "\t")
KM <- read.delim("PCSF/sscore_pcsf_prize_pval0.05_KM_diffWeights.txt", sep = "\t")

# WITH INTERACTOME
H9 <- openxlsx::read.xlsx("PCSF/sscorePLUSinteractome_pcsf_prize_pval0.05_diffWeights.xlsx", sheet = "H9")
JI <- openxlsx::read.xlsx("PCSF/sscorePLUSinteractome_pcsf_prize_pval0.05_diffWeights.xlsx", sheet = "JI")
KM <- openxlsx::read.xlsx("PCSF/sscorePLUSinteractome_pcsf_prize_pval0.05_diffWeights.xlsx", sheet = "KM")

# SET mu
mu = 0.0005

# CREATE PCSF NETWORKS
H9_pcsf_net <- makePCSFNetwork(H9, mu = mu)
JI_pcsf_net <- makePCSFNetwork(JI, mu = mu)
KM_pcsf_net <- makePCSFNetwork(KM, mu = mu)

# GENERATE PCSF PPI NETWORK
PCSFVisNodes(H9_pcsf_net, H9, name = paste0("H9_mu", mu))
PCSFVisNodes(JI_pcsf_net, JI, name = paste0("JI_mu", mu))
PCSFVisNodes(KM_pcsf_net, KM, name = paste0("KM_mu", mu))

# GENERATE PCSF INFLUENTIAL NODES NETWORK
PCSFVisInfluential(H9_pcsf_net, name = paste0("H9_mu", mu))
PCSFVisInfluential(JI_pcsf_net, name = paste0("JI_mu", mu))
PCSFVisInfluential(KM_pcsf_net, name = paste0("KM_mu", mu))

# GENERATE PCSF ENRICHED NETWORK
H9_pcsf_net_enriched <- pcsfRunEnrichment(H9_pcsf_net, name = paste0("H9_mu", mu))
pcsfEnrichedSubnet(H9_pcsf_net_enriched, name = paste0("H9_mu", mu))
JI_pcsf_net_enriched <- pcsfRunEnrichment(JI_pcsf_net, name = paste0("JI_mu", mu))
pcsfEnrichedSubnet(JI_pcsf_net_enriched, name = paste0("JI_mu", mu))
KM_pcsf_net_enriched <- pcsfRunEnrichment(KM_pcsf_net, name = paste0("KM_mu", mu))
pcsfEnrichedSubnet(KM_pcsf_net_enriched, name = paste0("KM_mu", mu))
