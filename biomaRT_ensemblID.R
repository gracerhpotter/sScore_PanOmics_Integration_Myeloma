
ensembl_dataset <- openxlsx::read.xlsx("sscoreInputs/RNAseq_LOGFC.xlsx")

ids <- ensembl_dataset %>% filter(is.na(X3))
ids <- ids$X1

mart <- biomaRt::useDataset("hsapiens_gene_ensembl", biomaRt::useMart("ensembl"))
ensembl <- biomaRt::useEnsembl(biomart = "ensembl", dataset = 'hsapiens_gene_ensembl')

gene_symbol_list <- biomaRt::getBM(filters= "ensembl_gene_id", 
                                   attributes= c("ensembl_gene_id", "external_gene_name", "uniprot_gn_symbol", "hgnc_symbol"), 
                                   values = ids, 
                                   mart = mart)

data_genes <- merge(ensembl_dataset, gene_symbol_list, by.x = "X1", by.y = "ensembl_gene_id") %>%
  filter(!(hgnc_symbol == "")) %>%
  dplyr::select(X1, hgnc_symbol)

data <- merge(ensembl_dataset, data_genes, by = "X1", all = TRUE)
openxlsx::write.xlsx(data, "sscoreInputs/RNAseq_LOGFC_biomaRt.xlsx")
