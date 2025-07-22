##############################################
# Differential Expression Analysis (DGEA) & GSEA
# Description:
#   - Perform DGEA using DESeq2 on filtered read counts
#   - Generate volcano plots and DEG tables
#   - Run GSEA with Hallmark gene sets using fgsea
##############################################

### Load Required Libraries ###
library(DESeq2)
library(ggplot2)
library(ggrepel)
library(openxlsx)

### Step 1: Prepare Input Count Matrix ###
# Define sample order
sample_order <- c("E01S1G", "E01S3G", "E01S5G",
                  "E03S1D", "E03S4D", "E04S4D", 
                  "E03S1G", "E03S4G", "E04S4G", 
                  "E06S1D", "E06S4D", "E06S3D", 
                  "E06S1G", "E06S4G", "E06S3G")

filtered_reads_matrix_counts <- reads_expr_filtered[, sample_order]

# Define sample conditions (one per column)
condition <- rep(c("Saline", "3x8Gy + Saline", 
                   "3x8Gy + Saline dist", 
                   "3x8Gy + Pembro", 
                   "3x8Gy + Pembro dist"), each = 3)

coldata <- data.frame(condition = condition)

# Round counts and create DESeq2 object
dds <- DESeqDataSetFromMatrix(
  countData = round(filtered_reads_matrix_counts),
  colData = coldata,
  design = ~ condition
)

# Run DESeq2
dds <- DESeq(dds)

### Step 2: Define DEG Processing Function ###
process_DEG <- function(deg_results, output_path, name) {
  deg_df <- as.data.frame(deg_results)
  
  # Classify DEGs
  deg_df$diffexpressed <- "NO"
  deg_df$diffexpressed[deg_df$log2FoldChange > 1 & deg_df$padj < 0.05] <- "UP"
  deg_df$diffexpressed[deg_df$log2FoldChange < -1 & deg_df$padj < 0.05] <- "DOWN"
  
  # Count significant DEGs
  cat("\n", name, "\n")
  cat("UP: ", sum(deg_df$diffexpressed == "UP"), "\n")
  cat("DOWN: ", sum(deg_df$diffexpressed == "DOWN"), "\n")
  
  # Highlight top 10 DEGs (5 up, 5 down)
  deg_df$gene_symbol <- rownames(deg_df)
  top_up <- head(deg_df[deg_df$diffexpressed == "UP", ][order(deg_df$padj), "gene_symbol"], 5)
  top_down <- head(deg_df[deg_df$diffexpressed == "DOWN", ][order(deg_df$padj), "gene_symbol"], 5)
  deg_df$delabel <- ifelse(deg_df$gene_symbol %in% c(top_up, top_down), deg_df$gene_symbol, NA)
  
  # Plot volcano
  p <- ggplot(deg_df, aes(x = log2FoldChange, y = -log10(padj), col = diffexpressed, label = delabel)) +
    geom_point() + 
    geom_text_repel(show.legend = FALSE) +
    scale_color_manual(values = c("DOWN" = "#009E73", "NO" = "black", "UP" = "#CC79A7")) +
    geom_vline(xintercept = c(-1, 1), col = "#CC79A7") +
    geom_hline(yintercept = -log10(0.05), col = "#CC79A7") +
    labs(title = name, y = "-log10(adj. p-value)", x = "log2 Fold Change") +
    theme_minimal(base_size = 15)
  
  # Save plot and results
  ggsave(file.path(output_path, paste0(name, "_volcano_plot.pdf")), plot = p, width = 6, height = 5, bg = "white")
  
  deg_filtered <- deg_df[deg_df$diffexpressed != "NO", ]
  deg_filtered <- deg_filtered[order(deg_filtered$log2FoldChange, decreasing = TRUE), ]
  write.xlsx(deg_filtered, file.path(output_path, paste0(name, "_DEG_list.xlsx")), row.names = FALSE)
}

### Step 3: Run All DEG Comparisons ###
deg_list <- list(
  RadioSaline_saline = results(dds, contrast = c("condition", "3x8Gy + Saline", "Saline")),
  RadioSalineDist_Saline = results(dds, contrast = c("condition", "3x8Gy + Saline dist", "Saline")),
  RadioSalineDist_RadioSaline = results(dds, contrast = c("condition", "3x8Gy + Saline dist", "3x8Gy + Saline")),
  RadioPembro_RadioSaline = results(dds, contrast = c("condition", "3x8Gy + Pembro", "3x8Gy + Saline")),
  RadioPembroDist_RadioSalineDist = results(dds, contrast = c("condition", "3x8Gy + Pembro dist", "3x8Gy + Saline dist")),
  RadioPembroDist_RadioPembro = results(dds, contrast = c("condition", "3x8Gy + Pembro dist", "3x8Gy + Pembro")),
  RadioPembro_Saline = results(dds, contrast = c("condition", "3x8Gy + Pembro", "Saline")),
  RadioPembroDist_Saline = results(dds, contrast = c("condition", "3x8Gy + Pembro dist", "Saline"))
)

output_path <- "results/dgea/"  # <-- Replace with your own relative output path

dir.create(output_path, recursive = TRUE, showWarnings = FALSE)
for (name in names(deg_list)) {
  process_DEG(deg_list[[name]], output_path, name)
}

### Step 4: GSEA Function ###
run_GSEA <- function(deg_list, pathways_file, output_path) {
  library(fgsea)
  library(tibble)
  library(dplyr)
  library(stringr)
  library(openxlsx)
  
  # Pathways to exclude
  remove_pathways <- c(
    "HALLMARK_PANCREAS_BETA_CELLS", "HALLMARK_SPERMATOGENESIS", 
    "HALLMARK_UV_RESPONSE_DN", "HALLMARK_UV_RESPONSE_UP", 
    "HALLMARK_ANDROGEN_RESPONSE", "HALLMARK_APICAL_SURFACE"
  )
  
  # Load Hallmark gene sets (GMT format)
  pathways <- gmtPathways(pathways_file)
  pathways <- pathways[!names(pathways) %in% remove_pathways]
  
  for (name in names(deg_list)) {
    message("Processing GSEA for: ", name)
    res <- as.data.frame(deg_list[[name]])
    res$padj <- res$padj + 1e-300  # Avoid log(0)
    res$ENTREZ <- rownames(res)
    res$metric <- -log10(res$padj) * sign(res$log2FoldChange)
    
    # Prepare ranking
    ranks <- deframe(res[, c("ENTREZ", "metric")])
    ranks <- na.omit(ranks)
    
    # Run fgsea
    set.seed(2)
    fgseaRes <- fgsea(pathways = pathways, stats = ranks)
    fgseaRes <- fgseaRes[order(fgseaRes$padj), ]
    
    # Save results
    write.xlsx(fgseaRes, file.path(output_path, paste0(name, "_GSEA_results.xlsx")))
    
    # Plot
    fgseaRes$signif <- ifelse(fgseaRes$padj < 0.25, "< 0.25", "> 0.25")
    fgseaRes$pathway <- gsub("_", " ", fgseaRes$pathway)
    fgseaRes$short <- ifelse(nchar(fgseaRes$pathway) > 40, paste0(substr(fgseaRes$pathway, 1, 39), "."), fgseaRes$pathway)
    
    p <- ggplot(fgseaRes, aes(reorder(short, NES), NES)) +
      geom_col(aes(fill = signif)) +
      scale_fill_manual(values = c("< 0.25" = "#009E73", "> 0.25" = "#999999")) +
      coord_flip() +
      labs(title = name, x = "Pathway", y = "NES") +
      theme_minimal()
    
    ggsave(file.path(output_path, paste0(name, "_GSEA_plot.pdf")), plot = p, width = 8, height = 6, bg = "white")
  }
}

### Step 5: Run GSEA ###
gsea_output <- "results/gsea/"
dir.create(gsea_output, recursive = TRUE, showWarnings = FALSE)

hallmark_pathways_file <- "data/pathways/h.all.v7.5.1.symbols.gmt.txt"  # <-- Replace with your general path
run_GSEA(deg_list, pathways_file = hallmark_pathways_file, output_path = gsea_output)
