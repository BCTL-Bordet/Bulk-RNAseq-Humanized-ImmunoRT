# This script processes human transcript quantification files,
# maps transcript IDs to gene IDs, aggregates expression data,
# and filters genes based on expression levels.

# Load required libraries
library(biomaRt)
library(dplyr)

#### Obtain human count matrices ####

# Define directories for input and output files
human_dir <- "folder_to_.sf_files"
output_dir <- "output_dir"

# Load Ensembl data to retrieve gene annotations
human_ensembl <- biomaRt::useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl", version = 112)

# Function to retrieve transcript-to-gene mappings
get_tx2gene <- function(ensembl) {
  getBM(attributes = c("ensembl_transcript_id", "ensembl_gene_id"), mart = ensembl)
}

# Retrieve the transcript-to-gene mapping for human
human_tx2gene <- get_tx2gene(human_ensembl)

# List all human sample files (adjust pattern if needed)
sample_files <- list.files(human_dir, pattern = "human_quant\\.sf$", full.names = TRUE)

# Initialize empty matrices for TPM and reads
tpm_matrix <- NULL
reads_matrix <- NULL

# Process each sample file
for (file in sample_files) {
  
  # Read the quantification file
  data <- read.table(file, header = FALSE, stringsAsFactors = FALSE)
  colnames(data) <- c("Name", "Length", "EffectiveLength", "TPM", "NumReads")
  
  # Remove version numbers from transcript IDs (e.g., ENST00000456328.2 -> ENST00000456328)
  data$Name <- sub("\\..*", "", data$Name)
  
  # Map transcript IDs to gene IDs
  gene_ids <- human_tx2gene[match(data$Name, human_tx2gene$ensembl_transcript_id), "ensembl_gene_id"]
  
  # Check for unmatched transcripts and print them
  if (any(is.na(gene_ids))) {
    unmatched_transcripts <- data$Name[is.na(gene_ids)]
    print("Unmatched transcripts:")
    print(unmatched_transcripts)
  }
  
  # Add the GeneID column to the data
  data$GeneID <- gene_ids
  
  # Aggregate TPM and NumReads by GeneID (sum values per gene)
  gene_data_tpm <- aggregate(data$TPM, by = list(data$GeneID), FUN = sum, na.rm = TRUE)
  gene_data_reads <- aggregate(data$NumReads, by = list(data$GeneID), FUN = sum, na.rm = TRUE)
  
  # Rename columns to include sample names
  colnames(gene_data_tpm) <- c("GeneID", sub(".*/(.*)_human_quant\\.sf", "\\1", file))
  colnames(gene_data_reads) <- c("GeneID", sub(".*/(.*)_human_quant\\.sf", "\\1", file))
  
  # Merge data with existing matrices
  if (is.null(tpm_matrix)) {
    tpm_matrix <- gene_data_tpm
    reads_matrix <- gene_data_reads
  } else {
    tpm_matrix <- merge(tpm_matrix, gene_data_tpm, by = "GeneID", all = TRUE)
    reads_matrix <- merge(reads_matrix, gene_data_reads, by = "GeneID", all = TRUE)
  }
}

# Save the aggregated matrices to CSV files
write.table(tpm_matrix, file = file.path(output_dir, "TPM_matrix.csv"), sep = ",", quote = FALSE, row.names = FALSE)
write.table(reads_matrix, file = file.path(output_dir, "Reads_matrix.csv"), sep = ",", quote = FALSE, row.names = FALSE)

#### Retrieve gene symbols and filter genes ####

# Retrieve gene symbols from Ensembl
gene_symbols_human <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"),
                            filters = "ensembl_gene_id",
                            values = reads_matrix$GeneID,
                            mart = human_ensembl)

# Map Gene Symbols to the matrices
reads_matrix$GeneSymbol <- gene_symbols_human$hgnc_symbol[match(reads_matrix$GeneID, gene_symbols_human$ensembl_gene_id)]
tpm_matrix$GeneSymbol <- gene_symbols_human$hgnc_symbol[match(tpm_matrix$GeneID, gene_symbols_human$ensembl_gene_id)]

# Remove rows with missing GeneSymbols
reads_matrix <- reads_matrix[!is.na(reads_matrix$GeneSymbol) & reads_matrix$GeneSymbol != "", ]

# Compute standard deviation and filter duplicated genes, keeping the highest variability
reads_matrix <- reads_matrix %>%
  group_by(GeneSymbol) %>%
  arrange(desc(apply(select(., -c(GeneID, GeneSymbol)), 1, sd))) %>%
  slice(1) %>%
  ungroup()

# Convert to data frame
reads_matrix <- as.data.frame(reads_matrix)

# Set row names based on GeneSymbol and remove redundant columns
rownames(reads_matrix) <- reads_matrix$GeneSymbol
reads_matrix$GeneID <- NULL
reads_matrix$GeneSymbol <- NULL

# Filter the TPM matrix to match filtered reads_matrix
tpm_matrix_filtered <- tpm_matrix[tpm_matrix$GeneID %in% rownames(reads_matrix), ]
rownames(tpm_matrix_filtered) <- tpm_matrix_filtered$GeneSymbol
tpm_matrix_filtered$GeneID <- NULL
tpm_matrix_filtered$GeneSymbol <- NULL

# Filter out lowly expressed genes (expression in at least 3 samples with count >= 10)
filtered_reads_matrix_counts <- reads_matrix[rowSums(reads_matrix >= 10) >= 3, ]
filtered_tpm_matrix_counts <- tpm_matrix_filtered[rownames(tpm_matrix_filtered) %in% rownames(filtered_reads_matrix_counts), ]

# Apply log transformation to TPM values
log_transformed_tpm <- log2(filtered_tpm_matrix_counts + 1)

