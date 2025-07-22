##############################################
# Pathway-Specific Heatmaps from Excel Sheets
# Description:
#   - Load gene sets from Excel (one pathway per sheet)
#   - Center gene expression by row (mean-centering only)
#   - Plot annotated heatmaps with treatment groups
##############################################

### Load Required Libraries ###
library(readxl)
library(randomcoloR)
library(ComplexHeatmap)

### Function to Process One Sheet ###
process_sheet_centering <- function(sheet_name, expression_matrix, file_path) {
  data <- read_excel(file_path, sheet = sheet_name)
  
  pathway_name <- unique(data$`Gene Set Name`)[1]
  genes <- data$Gene
  
  # Optional: annotation by gene if provided
  annotations <- if ("Annotation" %in% colnames(data)) {
    setNames(data$Annotation, data$Gene)
  } else {
    NULL
  }
  
  # Extract matching expression values
  selected_expression <- expression_matrix[match(genes, rownames(expression_matrix), nomatch = 0), ]
  selected_expression <- selected_expression[!is.na(rownames(selected_expression)), , drop = FALSE]
  
  if (nrow(selected_expression) == 0) {
    warning(paste("No matching genes found for pathway:", pathway_name))
    return(NULL)
  }
  
  # Mean-center expression values (z-score without scaling)
  rescaled_expression <- t(apply(selected_expression, 1, scale, center = TRUE, scale = FALSE))
  rownames(rescaled_expression) <- rownames(selected_expression)
  colnames(rescaled_expression) <- colnames(selected_expression)
  
  # Sample order
  column_order <- c(
    "E01S1G", "E01S3G", "E01S5G", 
    "E03S1D", "E03S4D", "E04S4D", 
    "E03S1G", "E03S4G", "E04S4G", 
    "E06S1D", "E06S4D", "E06S3D", 
    "E06S1G", "E06S4G", "E06S3G"
  )
  rescaled_expression <- rescaled_expression[, column_order, drop = FALSE]
  rescaled_expression <- rescaled_expression[genes[genes %in% rownames(rescaled_expression)], , drop = FALSE]
  
  # Column annotation (Treatment groups)
  treatment_levels <- c("Saline", "3x8Gy + Saline", "3x8Gy + Saline dist", "3x8Gy + Pembro", "3x8Gy + Pembro dist")
  treatment_colors <- c("Saline" = "lightblue", 
                        "3x8Gy + Saline" = "yellow3", 
                        "3x8Gy + Saline dist" = "purple4", 
                        "3x8Gy + Pembro" = "green4", 
                        "3x8Gy + Pembro dist" = "orange3")
  
  annotation_data <- data.frame(
    `Treatment group` = rep(treatment_levels, each = 3),
    row.names = column_order
  )
  
  column_annotation <- HeatmapAnnotation(
    "Treatment group" = annotation_data$`Treatment group`,
    col = list("Treatment group" = treatment_colors),
    annotation_legend_param = list(
      "Treatment group" = list(at = treatment_levels, title = "Treatment group")
    )
  )
  
  # Optional row annotation
  row_annotation <- if (!is.null(annotations)) {
    annotation_order <- unique(annotations)
    annotation_colors <- setNames(distinctColorPalette(length(annotation_order)), annotation_order)
    rowAnnotation(
      Annotation = annotations[rownames(rescaled_expression)],
      col = list(Annotation = annotation_colors),
      annotation_legend_param = list(
        Annotation = list(at = annotation_order, title = "Annotation")
      )
    )
  } else {
    NULL
  }
  
  # Plot the heatmap
  print(
    Heatmap(
      rescaled_expression,
      name = "log2(TPM)\nmean-centered",
      cluster_columns = FALSE,
      cluster_rows = FALSE,
      row_order = rownames(rescaled_expression),
      column_title = paste("Expression Heatmap -", pathway_name),
      row_title = "Genes",
      column_title_gp = gpar(fontsize = 12),
      row_title_gp = gpar(fontsize = 10),
      top_annotation = column_annotation,
      left_annotation = row_annotation
    )
  )
}

### Process All Sheets in the Excel File ###
process_all_sheets_centering <- function(file_path, expression_matrix) {
  sheet_names <- excel_sheets(file_path)
  for (sheet in sheet_names) {
    message("Processing sheet: ", sheet)
    process_sheet_centering(sheet, expression_matrix, file_path)
  }
}

### Example Usage ###
excel_path <- "data/pathways/gene_sets_by_pathway.xlsx"  # Generic relative path
process_all_sheets_centering(excel_path, log_transformed_tpm)
