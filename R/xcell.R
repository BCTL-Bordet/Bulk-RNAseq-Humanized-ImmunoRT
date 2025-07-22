##############################################
# xCell Deconvolution and Visualization
# Description:
#   - Perform cell-type deconvolution using xCell
#   - Filter out non-immune and irrelevant cell types
#   - Generate heatmap of rescaled xCell scores with treatment annotation
##############################################

### Load Required Libraries ###
library(xCell)
library(ComplexHeatmap)

### Step 1: Run xCell on log-transformed TPM ###
# Input: log2(TPM + 1) matrix (genes as rows, samples as columns)
xcell_scores <- xCellAnalysis(log_transformed_tpm)

# Inspect available cell types
all_cell_types <- rownames(xcell_scores)

# Define irrelevant or tissue-specific cell types to exclude
excluded_types <- c(
  "Astrocytes", "Chondrocytes", "Hepatocytes", "Keratinocytes", 
  "Melanocytes", "Mesangial cells", "Neurons", "Osteoblast",
  "Myocytes", "Skeletal muscle", "Smooth muscle", "Sebocytes",
  "ImmuneScore", "StromaScore", "MicroenvironmentScore", 
  "ly Endothelial cells", "Endothelial cells", "mv Endothelial cells", 
  "Adipocytes", "Epithelial cells", "Preadipocytes", "Pericytes", 
  "MSC", "Fibroblasts"
)

# Select only relevant immune cell types
cell_types_to_use <- setdiff(all_cell_types, excluded_types)

# Run xCell again using only the selected immune cell types
xcell_filtered <- as.data.frame(
  as.matrix(
    xCellAnalysis(log_transformed_tpm, cell.types.use = cell_types_to_use)
  )
)

### Step 2: Normalize Scores ###
# Row-wise z-score normalization
scaled_scores <- t(scale(t(xcell_filtered)))

### Step 3: Reorder Samples ###
# Define consistent sample order (based on your experimental design)
sample_order <- c(
  "E01S1G", "E01S3G", "E01S5G", 
  "E03S1D", "E03S4D", "E04S4D", 
  "E03S1G", "E03S4G", "E04S4G", 
  "E06S1D", "E06S4D", "E06S3D", 
  "E06S1G", "E06S4G", "E06S3G"
)

scaled_scores <- scaled_scores[, sample_order]

### Step 4: Prepare Column Annotation ###
# annotation_data must contain a column named 'Treatment group' (factors in proper order)
# Example:
# annotation_data <- data.frame(
#   `Treatment group` = rep(c("Saline", "3x8Gy + Saline", 
#                             "3x8Gy + Saline dist", 
#                             "3x8Gy + Pembro", 
#                             "3x8Gy + Pembro dist"), each = 3),
#   row.names = sample_order
# )

# Ensure correct order of factor levels
annotation_data$`Treatment group` <- factor(
  annotation_data$`Treatment group`,
  levels = c("Saline", "3x8Gy + Saline", 
             "3x8Gy + Saline dist", 
             "3x8Gy + Pembro", 
             "3x8Gy + Pembro dist")
)

# Define colors for annotation
column_annotation <- HeatmapAnnotation(
  "Treatment group" = annotation_data$`Treatment group`,
  col = list("Treatment group" = c(
    "Saline" = "lightblue",
    "3x8Gy + Saline" = "yellow3",
    "3x8Gy + Saline dist" = "purple4",
    "3x8Gy + Pembro" = "green4",
    "3x8Gy + Pembro dist" = "orange3"
  ))
)

### Step 5: Plot Heatmap ###
Heatmap(
  scaled_scores,
  top_annotation = column_annotation,
  name = "xCell deconvolution\n(rescaled)",
  cluster_columns = FALSE
)
