# This script processes human transcript quantification files,
# maps transcript IDs to gene IDs, aggregates expression data,
# and filters genes based on expression levels.

### Load Required Libraries ###
library(biomaRt)
library(dplyr)

### Define Input and Output Paths ###
human_dir <- "path/to/sf/files/"
output_dir <- "path/to/output/"

### Step 1: Retrieve Transcript-to-Gene Mapping from Ensembl ###
# Use specific Ensembl release for reproducibility
human_ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl", version = 112)

get_tx2gene <- function(ensembl) {
  getBM(attributes = c("ensembl_transcript_id", "ensembl_gene_id"), mart = ensembl)
}
human_tx2gene <- get_tx2gene(human_ensembl)

### Step 2: Load Salmon Quantifications and Aggregate to Gene Level ###
sample_files <- list.files(human_dir, pattern = "human_quant\\.sf$", full.names = TRUE)
tpm_matrix <- NULL
reads_matrix <- NULL

for (file in sample_files) {
  data <- read.table(file, header = FALSE, stringsAsFactors = FALSE)
  colnames(data) <- c("Name", "Length", "EffectiveLength", "TPM", "NumReads")
  
  # Remove version numbers from transcript IDs (e.g., ENST00000456328.2 -> ENST00000456328)
  data$Name <- sub("\\..*", "", data$Name)
  
  # Map transcripts to genes
  data$GeneID <- human_tx2gene$ensembl_gene_id[match(data$Name, human_tx2gene$ensembl_transcript_id)]
  
  # Print unmatched transcript IDs (for debugging)
  if (any(is.na(data$GeneID))) {
    cat("Unmatched transcripts:\n")
    print(data$Name[is.na(data$GeneID)])
  }
  
  # Aggregate TPM and reads at the gene level
  sample_name <- sub(".*/(.*)_human_quant\\.sf", "\\1", file)
  gene_data_tpm <- aggregate(TPM ~ GeneID, data, sum, na.rm = TRUE)
  gene_data_reads <- aggregate(NumReads ~ GeneID, data, sum, na.rm = TRUE)
  
  colnames(gene_data_tpm) <- c("GeneID", sample_name)
  colnames(gene_data_reads) <- c("GeneID", sample_name)
  
  # Merge with main matrices
  tpm_matrix <- if (is.null(tpm_matrix)) gene_data_tpm else merge(tpm_matrix, gene_data_tpm, by = "GeneID", all = TRUE)
  reads_matrix <- if (is.null(reads_matrix)) gene_data_reads else merge(reads_matrix, gene_data_reads, by = "GeneID", all = TRUE)
}

# Export raw matrices
write.table(tpm_matrix, file = file.path(output_dir, "TPM_matrix.csv"), sep = ",", quote = FALSE, row.names = FALSE)
write.table(reads_matrix, file = file.path(output_dir, "Reads_matrix.csv"), sep = ",", quote = FALSE, row.names = FALSE)

### Step 3: Annotate Gene Symbols ###
# Reload Ensembl for updated symbols
human_ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

gene_symbols_human <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"),
                            filters = "ensembl_gene_id",
                            values = reads_matrix$GeneID,
                            mart = human_ensembl)

# Map gene symbols to matrices
reads_matrix$GeneSymbol <- gene_symbols_human$hgnc_symbol[match(reads_matrix$GeneID, gene_symbols_human$ensembl_gene_id)]
tpm_matrix$GeneSymbol <- gene_symbols_human$hgnc_symbol[match(tpm_matrix$GeneID, gene_symbols_human$ensembl_gene_id)]

# Filter out genes with no GeneSymbol
tpm_matrix <- tpm_matrix[!is.na(tpm_matrix$GeneSymbol) & tpm_matrix$GeneSymbol != "", ]
reads_matrix <- reads_matrix[!is.na(reads_matrix$GeneSymbol) & reads_matrix$GeneSymbol != "", ]

### Step 4: Deduplicate by Gene Symbol Using Standard Deviation ###
# Compute SD across TPM values
tpm_numeric <- tpm_matrix[, !(names(tpm_matrix) %in% c("GeneID", "GeneSymbol"))]
tpm_matrix$std_dev <- apply(tpm_numeric, 1, sd)

# Keep the entry with highest SD for each GeneSymbol
tpm_matrix_dedup <- tpm_matrix %>%
  group_by(GeneSymbol) %>%
  arrange(desc(std_dev)) %>%
  slice(1) %>%
  ungroup() %>%
  as.data.frame()

# Filter reads matrix accordingly
gene_ids <- tpm_matrix_dedup$GeneID
reads_matrix_filtered <- reads_matrix[reads_matrix$GeneID %in% gene_ids, ]
tpm_matrix_filtered   <- tpm_matrix_dedup

### Step 5: Filter Lowly Expressed Genes ###
# Remove genes with <10 reads in fewer than 3 samples
reads_expr <- reads_matrix_filtered[, !(names(reads_matrix_filtered) %in% c("GeneID", "GeneSymbol"))]
rownames(reads_expr) <- reads_matrix_filtered$GeneSymbol

reads_expr_filtered <- reads_expr[rowSums(reads_expr >= 10) >= 3, ]

### Step 6: Align TPM Matrix and Log-Transform ###
tpm_expr <- tpm_matrix_filtered[, !(names(tpm_matrix_filtered) %in% c("GeneID", "GeneSymbol", "std_dev"))]
rownames(tpm_expr) <- tpm_matrix_filtered$GeneSymbol

tpm_expr_filtered <- tpm_expr[rownames(tpm_expr) %in% rownames(reads_expr_filtered), ]
log_transformed_tpm <- log2(tpm_expr_filtered + 1)

### Output ###
# Final log-transformed TPM matrix (genes x samples) is stored in:
#   log_transformed_tpm

# You may save it for downstream use:
write.table(log_transformed_tpm, file = file.path(output_dir, "log2_TPM_matrix.csv"),
            sep = ",", quote = FALSE, row.names = TRUE)
  
