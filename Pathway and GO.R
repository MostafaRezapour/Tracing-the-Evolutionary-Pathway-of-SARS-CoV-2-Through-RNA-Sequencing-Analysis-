print("R Script is running ...")
cat("\014")
rm(list = ls())
getwd()
setwd(getwd())
library(enrichR)
library(clusterProfiler)
library(org.Hs.eg.db)
library(AnnotationDbi)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("No input file provided. Please provide a CSV file path as the first argument.")
}
csv_file <- args[1]  # CSV file path
output_base <- sub("\\.csv$", "", basename(csv_file))  # Extract base name without .csv
if (!file.exists(csv_file)) {
  stop(paste("The file", csv_file, "does not exist."))
}
genes_data <- read.csv(csv_file, stringsAsFactors = FALSE)
upregulated_genes <- genes_data$Gene.Symbol[genes_data$Regulation == "Upregulated"]
downregulated_genes <- genes_data$Gene.Symbol[genes_data$Regulation == "Downregulated"]
upregulated_genes <- unique(na.omit(upregulated_genes))
downregulated_genes <- unique(na.omit(downregulated_genes))
map_genes_to_entrez <- function(genes) {
  mapped <- bitr(genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
  mapped <- mapped[!is.na(mapped$ENTREZID), ]  # Remove unmapped genes
  mapped <- mapped[!duplicated(mapped$SYMBOL), ]  # Remove duplicates
  return(mapped$ENTREZID)
}
upregulated_entrez <- map_genes_to_entrez(upregulated_genes)
downregulated_entrez <- map_genes_to_entrez(downregulated_genes)
perform_go_analysis <- function(entrez_ids, output_file) {
  cat("Performing GO analysis for:", output_file, "\n")
  cat("Number of ENTREZ IDs:", length(entrez_ids), "\n")
  if (length(entrez_ids) == 0) {
    cat("No genes available for GO enrichment analysis in", output_file, "\n")
    return(NULL)
  }
 
  go_results <- enrichGO(
    gene = entrez_ids,
    OrgDb = org.Hs.eg.db,
    ont = "BP",  # Biological Processes
    pAdjustMethod = "BH",  # Benjamini-Hochberg correction
    qvalueCutoff = 0.05
  )
  
  if (!is.null(go_results) && nrow(as.data.frame(go_results)) > 0) {
    cat("Significant GO terms identified for", output_file, "\n")
    go_df <- as.data.frame(go_results)
    go_df$GeneSymbols <- sapply(go_df$geneID, function(ids) {
      entrez_ids <- strsplit(ids, "/")[[1]]
      symbols <- bitr(entrez_ids, fromType = "ENTREZID", toType = "SYMBOL", OrgDb = org.Hs.eg.db)
      paste(symbols$SYMBOL, collapse = ", ")
    })
    write.csv(go_df, output_file, row.names = FALSE)
    cat("GO enrichment results saved to", output_file, "\n")
  } else {
    cat("No significant GO terms found for", output_file, "\n")
  }
}
perform_go_analysis(upregulated_entrez, paste0("GO_Upregulated_", output_base, ".csv"))
perform_go_analysis(downregulated_entrez, paste0("GO_Downregulated_", output_base, ".csv"))
cat("R Script execution completed successfully.\n")




packages <- c("clusterProfiler", "org.Hs.eg.db", "AnnotationDbi", "enrichR")
for (pkg in packages) {
  if (requireNamespace(pkg, quietly = TRUE)) {
    cat(pkg, "version:", as.character(packageVersion(pkg)), "\n")
  } else {
    cat(pkg, "is not installed.\n")
  }
}

