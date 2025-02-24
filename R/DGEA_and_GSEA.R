#### Differential Gene Expression Analysis (DGEA) ####

# Load necessary libraries
library(DESeq2)
library(pheatmap)
library(dendextend)
library(ggplot2)
library(ggrepel)
library(openxlsx)

# Reorder columns of the count matrix based on sample names
column_order <- c("E01S1G", "E01S3G", "E01S5G", "E02S2G", "E02S3G", "E02S4G", 
                  "E03S1D", "E03S4D", "E04S4D", "E03S1G", "E03S4G", "E04S4G", 
                  "E06S1D", "E06S4D", "E06S3D", "E06S1G", "E06S4G", "E06S3G")
filtered_reads_matrix_counts <- filtered_reads_matrix_counts[, column_order]

# Define experimental conditions for each sample
treatment_conditions <- c("Saline", "Pembro", "3x8Gy + Saline", 
                          "3x8Gy + Saline dist", "3x8Gy + Pembro", "3x8Gy + Pembro dist")
condition <- rep(treatment_conditions, each = 3)
coldata <- data.frame(condition = condition)

# Create DESeq2 dataset
dds <- DESeqDataSetFromMatrix(
  countData = round(filtered_reads_matrix_counts),
  colData = coldata,
  design = ~ condition
)

# Normalize counts using variance stabilizing transformation (VST)
dds <- estimateSizeFactors(dds)
normalized_counts <- log2(counts(dds, normalized = TRUE) + 1e-15)

#### Data Visualization ####

# Perform variance stabilizing transformation
vsd <- vst(dds, blind = FALSE)

# Compute sample distance matrix
dist_matrix <- dist(t(assay(vsd)))
dist_matrix <- as.matrix(dist_matrix)
dist_obj <- as.dist(dist_matrix)

# Perform hierarchical clustering
hc <- hclust(dist_obj)

# Define colors for different conditions
class_colors <- c("Saline" = "lightblue", "Pembro" = "pink", "3x8Gy + Saline" = "yellow3", 
                  "3x8Gy + Saline dist" = "purple4", "3x8Gy + Pembro" = "green4", "3x8Gy + Pembro dist" = "orange3")
names(class_colors) <- treatment_conditions

# Convert hierarchical clustering to dendrogram
dend <- as.dendrogram(hc) %>%
  set("labels_colors", value = class_colors[condition[labels(dend)]]) %>%
  set("labels_cex", value = 0.8)

# Plot hierarchical clustering dendrogram
par(mar = c(5, 5, 4, 10), xpd = TRUE)
plot(dend, main = "Hierarchical Clustering Dendrogram",
     xlab = "Samples", ylab = "Distance")
legend("topright", legend = names(class_colors), col = class_colors, pch = 16, 
       title = "Treatment group", bty = "n", cex = 0.7, x.intersp = 0.6, y.intersp = 0.6, inset = c(-0.55, -0.2))

# Generate heatmap of sample distances
pheatmap(as.matrix(dist_matrix), clustering_distance_rows = dist_matrix, clustering_distance_cols = dist_matrix)

# Perform PCA analysis
pca_data <- plotPCA(vsd, intgroup = "condition", returnData = TRUE)
pca_data$condition <- factor(pca_data$condition, levels = treatment_conditions)

# Plot PCA
ggplot(pca_data, aes(x = PC1, y = PC2, color = condition, label = name)) +
  geom_point(size = 4, alpha = 0.8) +
  scale_color_manual(values = class_colors) +
  theme_minimal() +
  labs(title = "PCA Plot", x = "PC1", y = "PC2", color = "Treatment group") +
  theme(legend.title = element_text(face = "bold"))

#### Differential Expression Analysis ####

dds <- DESeq(dds)
results <- results(dds)

# Function to generate volcano plots and save differentially expressed genes (DEGs)
process_DEG <- function(deg_results, name) {
  deg_df <- as.data.frame(deg_results)
  
  # Define differentially expressed genes
  deg_df$diffexpressed <- "NO"
  deg_df$diffexpressed[deg_df$log2FoldChange > 1 & deg_df$padj < 0.05] <- "UP"
  deg_df$diffexpressed[deg_df$log2FoldChange < -1 & deg_df$padj < 0.05] <- "DOWN"
  
  # Count significant DEGs
  num_up <- sum(deg_df$diffexpressed == "UP")
  num_down <- sum(deg_df$diffexpressed == "DOWN")
  
  print(paste(name, "- Upregulated genes:", num_up, "| Downregulated genes:", num_down))
  
  # Select top significant genes
  deg_df$gene_symbol <- rownames(deg_df)
  top_up <- head(deg_df[deg_df$diffexpressed == "UP", ][order(deg_df$padj[deg_df$diffexpressed == "UP"], na.last = NA), "gene_symbol"], 5)
  top_down <- head(deg_df[deg_df$diffexpressed == "DOWN", ][order(deg_df$padj[deg_df$diffexpressed == "DOWN"], na.last = NA), "gene_symbol"], 5)
  deg_df$delabel <- ifelse(deg_df$gene_symbol %in% c(top_up, top_down), deg_df$gene_symbol, NA)
  
  # Create volcano plot
  p <- ggplot(data = deg_df, aes(x = log2FoldChange, y = -log10(padj), col = diffexpressed, label = delabel)) +
    geom_point() +
    theme_minimal() +
    geom_text_repel(show.legend = FALSE) +
    scale_color_manual(values = c("#CC79A7", "black", "#009E73")) +
    geom_vline(xintercept = c(-1, 1), col = "#CC79A7") +
    geom_hline(yintercept = -log10(0.05), col = "#CC79A7") +
    labs(y = "-log10(Adj. P value)", x = "log2 Fold Change") +
    theme(axis.text = element_text(color = "black"), text = element_text(size = 15))
  
  print(p)
}

# Define contrasts for DEG analysis
deg_list <- list(
  pembro_saline = results(dds, contrast = c("condition", "Pembro", "Saline")),
  RadioSaline_saline = results(dds, contrast = c("condition", "3x8Gy + Saline", "Saline")),
  RadioSalineDist_Saline = results(dds, contrast = c("condition", "3x8Gy + Saline dist", "Saline")),
  PembroRadio_pembro = results(dds, contrast = c("condition", "3x8Gy + Pembro", "Pembro")),
  PembroRadioDist_pembro = results(dds, contrast = c("condition", "3x8Gy + Pembro dist", "Pembro"))
)

# Process each DEG comparison
for (name in names(deg_list)) {
  process_DEG(deg_list[[name]], name)
}
#### Gene Set Enrichment Analysis (GSEA) ####
run_GSEA <- function(deg_list, pathways_file, output_path) {
  library(fgsea)
  library(tibble)
  library(ggplot2)
  library(openxlsx)
  library(dplyr)
  library(stringr)
  
  # Define pathways to exclude from analysis
  remove_pathways <- c(
    "HALLMARK_PANCREAS_BETA_CELLS", "HALLMARK_SPERMATOGENESIS", 
    "HALLMARK_UV_RESPONSE_DN", "HALLMARK_UV_RESPONSE_UP", 
    "HALLMARK_ANDROGEN_RESPONSE", "HALLMARK_APICAL_SURFACE"
  )
  
  # Load gene set pathways
  pathways.hallmark <- gmtPathways(pathways_file)
  
  # Remove unwanted pathways
  pathways.hallmark <- pathways.hallmark[!names(pathways.hallmark) %in% remove_pathways]
  
  # Loop through each DEG dataset
  for (name in names(deg_list)) {
    message("Processing: ", name)
    
    # Extract DEG results
    deg_results <- deg_list[[name]]
    
    # Compute ranking metric
    deg_results$p_val2 <- deg_results$padj + 1e-300  # Prevent log(0)
    deg_results$ENTREZ <- rownames(deg_results)
    deg_results$fcsign <- sign(deg_results$log2FoldChange)
    deg_results$logP <- -log10(deg_results$p_val2)
    deg_results$metric <- deg_results$logP / deg_results$fcsign
    
    # Prepare ranking table
    res2 <- deg_results[, c("ENTREZ", "metric")]
    colnames(res2) <- c("SYMBOL", "metric")
    res2 <- na.omit(res2)
    ranks <- deframe(res2)
    
    # Perform Gene Set Enrichment Analysis (GSEA)
    set.seed(2)
    fgseaRes <- fgsea(pathways = pathways.hallmark, stats = ranks)
    
    # Sort results by adjusted p-value
    fgseaResTidy <- fgseaRes[order(fgseaRes$padj), ]
    
    # Save results to an Excel file
    excel_file <- paste0(output_path, name, "_GSEA_results.xlsx")
    write.xlsx(fgseaResTidy, excel_file)
    
    # Mark significantly enriched pathways
    fgseaResTidy$`p.adjust` <- ifelse(fgseaResTidy$padj < 0.25, "< 0.25", "> 0.25")
    group.colors <- c(`< 0.25` = "#009E73", `> 0.25` = "#999999")
    fgseaResTidy$pathway <- gsub("_", " ", fgseaResTidy$pathway)
    
    # Shorten long pathway names
    fgseaResTidy$pathwaysss <- sapply(fgseaResTidy$pathway, function(x) {
      if (str_length(x) > 40) paste0(substr(x, 1, 39), ".") else x
    })
    
    # Generate GSEA plot
    p <- ggplot(fgseaResTidy, aes(reorder(pathwaysss, NES), NES)) +
      geom_col(aes(fill = `p.adjust`)) + 
      scale_fill_manual(values = group.colors) +
      coord_flip() +
      labs(x = "Pathway", y = "NES", title = name) +
      theme(
        axis.text = element_text(colour = "black"), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black")
      )
    
    # Save plot as PDF
    pdf_file <- paste0(output_path, name, "_GSEA_plot.pdf")
    ggsave(pdf_file, plot = p, width = 8, height = 6, bg = "white")
    
    message("Saved: ", pdf_file, " and ", excel_file)
  }
}

# Define input files and directories for GSEA analysis
output_directory <- "your/output/directory/"
pathways.hallmark <- "your/pathways/file.gmt"

# Run GSEA analysis on DEG results
run_GSEA(
  deg_list = deg_list,  # List of differentially expressed genes
  pathways_file = pathways.hallmark,  # Path to hallmark gene sets
  output_path = output_directory  # Output directory for results
)
