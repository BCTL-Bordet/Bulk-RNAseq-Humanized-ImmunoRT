#### xCell analysis ####
# Load required libraries
library(xCell)
library(ComplexHeatmap)

# Perform xCell analysis on log-transformed TPM matrix
xcell_human <- xCellAnalysis(log_transformed_tpm)

# Extract cell type names
cell_types <- rownames(xcell_human)

# Define cell types to exclude from the analysis
cell_types_to_exclude <- c("Astrocytes", "Chondrocytes", "Hepatocytes", "Keratinocytes", 
                           "Melanocytes", "Mesangial cells", "Neurons", "Osteoblast", 
                           "Myocytes", "Skeletal muscle", "Smooth muscle", "Sebocytes", 
                           "ImmuneScore", "StromaScore", "MicroenvironmentScore", 
                           "ly Endothelial cells", "Endothelial cells", 
                           "mv Endothelial cells", "Adipocytes", "Epithelial cells", 
                           "Preadipocytes", "Pericytes", "MSC", "Fibroblasts")

# Filter cell types to include only relevant ones
cell_types_to_use <- cell_types[!cell_types %in% cell_types_to_exclude]

# Perform xCell analysis again with selected cell types
table_xcell <- as.data.frame(as.matrix(xCellAnalysis(log_transformed_tpm, cell.types.use = cell_types_to_use)))

# Scale scores across samples
scaled_scores <- t(scale(t(table_xcell)))  # Normalize per cell type

# Define specific order for sample names
sample_order <- c("E01S1G", "E01S3G", "E01S5G", "E02S2G", "E02S3G", "E02S4G", 
                  "E03S1D", "E03S4D", "E04S4D", "E03S1G", "E03S4G", "E04S4G", 
                  "E06S1D", "E06S4D", "E06S3D", "E06S1G", "E06S4G", "E06S3G")

# Reorder scaled scores by predefined sample order
scaled_scores <- scaled_scores[, sample_order]

# Define annotation for treatment groups
annotation_data <- data.frame(
  Group = rep(c("Saline", "Pembro", "3x8Gy + Saline", "3x8Gy + Saline dist", 
                "3x8Gy + Pembro", "3x8Gy + Pembro dist"), each = 3)
)
rownames(annotation_data) <- colnames(scaled_scores)
colnames(annotation_data) <- "Treatment group"

# Create heatmap annotation object
column_annotation <- HeatmapAnnotation(
  "Treatment group" = annotation_data$`Treatment group`,
  col = list(
    "Treatment group" = c("Saline" = "lightblue", "Pembro" = "pink", 
                          "3x8Gy + Saline" = "yellow3", "3x8Gy + Saline dist" = "purple4", 
                          "3x8Gy + Pembro" = "green4", "3x8Gy + Pembro dist" = "orange3")
  )
)

# Generate heatmap visualization
Heatmap(
  scaled_scores,
  top_annotation = column_annotation,
  name = "xCell deconvolution\nrescaled", 
  cluster_columns = FALSE
)

#### Gene expression heatmaps (complete) ####
# Load required libraries
library(readxl)
library(randomcoloR)
library(ComplexHeatmap)

# Function to process each sheet of an Excel file
process_sheet <- function(sheet_name, expression_matrix, file_path) {
  # Read the sheet containing pathway and gene information
  data <- read_excel(file_path, sheet = sheet_name)
  
  # Extract pathway name and associated genes
  pathway_name <- unique(data$`Gene Set Name`)[1]
  genes <- data$Gene
  
  # Check for the presence of "Annotation" column
  if ("Annotation" %in% colnames(data)) {
    annotations <- setNames(data$Annotation, data$Gene)  # Assign gene names as names
  } else {
    annotations <- NULL
  }
  
  # Subset the expression matrix based on selected genes while preserving the original order
  selected_expression <- expression_matrix[match(genes, rownames(expression_matrix), nomatch = 0), ]
  
  # Remove NA rows (genes not found in the expression matrix)
  selected_expression <- selected_expression[!is.na(rownames(selected_expression)), , drop = FALSE]
  
  # Check if there are any matched genes
  if (nrow(selected_expression) == 0) {
    warning(paste("No matching genes found in the expression matrix for pathway:", pathway_name))
    return(NULL)
  }
  
  # Normalize expression values using z-score transformation
  rescaled_expression <- t(apply(selected_expression, 1, scale))
  rownames(rescaled_expression) <- rownames(selected_expression)
  colnames(rescaled_expression) <- colnames(selected_expression)
  
  # Define the order of samples for visualization
  column_order <- c("E01S1G", "E01S3G", "E01S5G", "E02S2G", "E02S3G", "E02S4G", 
                    "E03S1D", "E03S4D", "E04S4D", "E03S1G", "E03S4G", "E04S4G", 
                    "E06S1D", "E06S4D", "E06S3D", "E06S1G", "E06S4G", "E06S3G")
  
  # Reorder the expression matrix columns
  rescaled_expression <- rescaled_expression[, column_order, drop = FALSE]
  
  # Preserve the original order of genes from the Excel sheet
  rescaled_expression <- rescaled_expression[genes[genes %in% rownames(rescaled_expression)], , drop = FALSE]
  
  # Create annotation for treatment groups
  annotation_data <- data.frame(
    Group = rep(c("Saline", "Pembro", "3x8Gy + Saline", "3x8Gy + Saline dist", 
                  "3x8Gy + Pembro", "3x8Gy + Pembro dist"), each = 3)
  )
  rownames(annotation_data) <- column_order
  colnames(annotation_data) <- "Treatment group"
  
  # Define heatmap column annotation
  column_annotation <- HeatmapAnnotation(
    "Treatment group" = annotation_data$`Treatment group`,
    col = list(
      "Treatment group" = c("Saline" = "lightblue", "Pembro" = "pink", 
                            "3x8Gy + Saline" = "yellow3", "3x8Gy + Saline dist" = "purple4", 
                            "3x8Gy + Pembro" = "green4", "3x8Gy + Pembro dist" = "orange3")
    )
  )
  
  # Define row annotation if "Annotation" column exists
  if (!is.null(annotations)) {
    unique_annotations <- unique(annotations)
    annotation_colors <- setNames(distinctColorPalette(length(unique_annotations)), unique_annotations)
    
    row_annotation <- rowAnnotation(
      Annotation = annotations[rownames(rescaled_expression)],
      col = list(Annotation = annotation_colors)
    )
  } else {
    row_annotation <- NULL
  }
  
  # Generate heatmap visualization
  print(Heatmap(
    rescaled_expression, 
    name = "Rescaled\nexpression",
    cluster_columns = FALSE,  # Disable column clustering
    cluster_rows = FALSE,  # Disable row clustering
    row_order = rownames(rescaled_expression),  # Maintain gene order
    column_title = paste("Expression Heatmap -", pathway_name),
    row_title = "Genes", 
    column_title_gp = gpar(fontsize = 12),
    row_title_gp = gpar(fontsize = 10),
    top_annotation = column_annotation,
    left_annotation = row_annotation
  ))
}

# Function to process all sheets in an Excel file
process_all_sheets <- function(file_path, expression_matrix) {
  # Retrieve sheet names
  sheet_names <- excel_sheets(file_path)
  
  # Iterate over each sheet and process it
  for (sheet in sheet_names) {
    print(paste("Processing sheet:", sheet))
    process_sheet(sheet, expression_matrix, file_path)
  }
}

# Define file path and execute processing
file_path <- "path_to_file/25_02_05_RNAseq_Final_Heatmaps_Clean_V2.xlsx"  # Replace with actual file path
process_all_sheets(file_path, log_transformed_tpm)



#### Gene expression heatmaps (excluding Saline and Pembro samples) ####
# Load required libraries
library(readxl)
library(randomcoloR)
library(ComplexHeatmap)

# Function to process each sheet of an Excel file
process_sheet <- function(sheet_name, expression_matrix, file_path) {
  # Read the sheet containing pathway and gene information
  data <- read_excel(file_path, sheet = sheet_name)
  
  # Extract pathway name and associated genes
  pathway_name <- unique(data$`Gene Set Name`)[1]
  genes <- data$Gene
  
  # Check for the presence of "Annotation" column
  if ("Annotation" %in% colnames(data)) {
    annotations <- setNames(data$Annotation, data$Gene)  # Assign gene names as names
  } else {
    annotations <- NULL
  }
  
  # Subset the expression matrix based on selected genes while preserving the original order
  selected_expression <- expression_matrix[match(genes, rownames(expression_matrix), nomatch = 0), ]
  
  # Remove NA rows (genes not found in the expression matrix)
  selected_expression <- selected_expression[!is.na(rownames(selected_expression)), , drop = FALSE]
  
  # Check if there are any matched genes
  if (nrow(selected_expression) == 0) {
    warning(paste("No matching genes found in the expression matrix for pathway:", pathway_name))
    return(NULL)
  }
  
  # Define order of all samples
  column_order <- c("E01S1G", "E01S3G", "E01S5G", "E02S2G", "E02S3G", "E02S4G", 
                    "E03S1D", "E03S4D", "E04S4D", "E03S1G", "E03S4G", "E04S4G", 
                    "E06S1D", "E06S4D", "E06S3D", "E06S1G", "E06S4G", "E06S3G")
  
  # Select a subset of relevant samples for visualization
  columns_to_include <- c("E03S1D", "E03S4D", "E04S4D", "E03S1G", "E03S4G", "E04S4G", 
                          "E06S1D", "E06S4D", "E06S3D", "E06S1G", "E06S4G", "E06S3G")
  
  # Reorder and filter selected expression matrix
  selected_expression <- selected_expression[, columns_to_include, drop = FALSE]
  
  # Normalize expression values using z-score transformation
  rescaled_expression <- t(apply(selected_expression, 1, scale))
  rownames(rescaled_expression) <- rownames(selected_expression)
  colnames(rescaled_expression) <- colnames(selected_expression)
  
  # Preserve the original order of genes from the Excel sheet
  rescaled_expression <- rescaled_expression[genes[genes %in% rownames(rescaled_expression)], , drop = FALSE]
  
  # Create annotation for treatment groups
  annotation_data <- data.frame(
    Group = rep(c("3x8Gy + Saline", "3x8Gy + Saline dist", "3x8Gy + Pembro", "3x8Gy + Pembro dist"), each = 3)
  )
  rownames(annotation_data) <- columns_to_include
  colnames(annotation_data) <- "Treatment group"
  
  # Define heatmap column annotation
  column_annotation <- HeatmapAnnotation(
    "Treatment group" = annotation_data$`Treatment group`,
    col = list(
      "Treatment group" = c("3x8Gy + Saline" = "yellow3", "3x8Gy + Saline dist" = "purple4", 
                            "3x8Gy + Pembro" = "green4", "3x8Gy + Pembro dist" = "orange3")
    )
  )
  
  # Define row annotation if "Annotation" column exists
  if (!is.null(annotations)) {
    unique_annotations <- unique(annotations)
    annotation_colors <- setNames(distinctColorPalette(length(unique_annotations)), unique_annotations)
    
    row_annotation <- rowAnnotation(
      Annotation = annotations[rownames(rescaled_expression)],
      col = list(Annotation = annotation_colors)
    )
  } else {
    row_annotation <- NULL
  }
  
  # Generate heatmap visualization
  print(Heatmap(
    rescaled_expression, 
    name = "Rescaled\nexpression",
    cluster_columns = FALSE,  # Disable column clustering
    cluster_rows = FALSE,  # Disable row clustering
    row_order = rownames(rescaled_expression),  # Maintain gene order
    column_title = paste("Expression Heatmap -", pathway_name),
    row_title = "Genes", 
    column_title_gp = gpar(fontsize = 12),
    row_title_gp = gpar(fontsize = 10),
    top_annotation = column_annotation,
    left_annotation = row_annotation
  ))
}

# Function to process all sheets in an Excel file
process_all_sheets <- function(file_path, expression_matrix) {
  # Retrieve sheet names
  sheet_names <- excel_sheets(file_path)
  
  # Iterate over each sheet and process it
  for (sheet in sheet_names) {
    print(paste("Processing sheet:", sheet))
    process_sheet(sheet, expression_matrix, file_path)
  }
}

# Define file path and execute processing
file_path <- "path_to_file/25_02_05_RNAseq_Final_Heatmaps_Clean_V2.xlsx"  # Replace with actual file path
process_all_sheets(file_path, log_transformed_tpm)
