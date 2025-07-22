##############################################
# GSVA Analysis
# Description:
#   - Perform GSVA using Hallmark gene sets from MSigDB
#   - Filter biologically irrelevant pathways
#   - Plot heatmap of GSVA scores with treatment group annotation
##############################################

### Load Required Libraries ###
library(GSVA)
library(dplyr)
library(ComplexHeatmap)

### Step 1: Load Gene Sets ###
# Load MSigDB gene sets (preprocessed as a data frame)
all_gene_sets <- readRDS("/path/to/file.gmt")

# Filter for Hallmark gene sets (category "H")
h_gene_sets <- all_gene_sets %>%
  filter(gs_cat == "H") %>%
  filter(!gs_name %in% c(
    "HALLMARK_PANCREAS_BETA_CELLS", 
    "HALLMARK_SPERMATOGENESIS", 
    "HALLMARK_UV_RESPONSE_DN", 
    "HALLMARK_UV_RESPONSE_UP", 
    "HALLMARK_ANDROGEN_RESPONSE", 
    "HALLMARK_APICAL_SURFACE"
  ))

# Display the retained gene set names
unique(h_gene_sets$gs_name)

# Prepare list of gene sets for GSVA (gene symbols per pathway)
msigdbr_list <- split(h_gene_sets$gene_symbol, f = h_gene_sets$gs_name)
msigdbr_list_unique <- lapply(msigdbr_list, unique)

### Step 2: Run GSVA ###
# Input: log2-transformed TPM matrix (genes as rows, samples as columns)
params <- gsvaParam(
  as.matrix(log_transformed_tpm), 
  geneSets = msigdbr_list_unique, 
  kcdf = "Gaussian"
)

# Compute GSVA scores and transpose (samples as rows)
gsva_res_filtered <- as.data.frame(t(gsva(params)))

# Clean up gene set names (remove "HALLMARK_" and underscores)
colnames(gsva_res_filtered) <- colnames(gsva_res_filtered) %>%
  gsub("HALLMARK_", "", .) %>%
  gsub("_", " ", .)

### Step 3: Subset and Reorder Samples (if needed) ###
# Define specific sample order for consistent plotting
sample_order <- c(
  "E01S1G", "E01S3G", "E01S5G", 
  "E03S1D", "E03S4D", "E04S4D", 
  "E03S1G", "E03S4G", "E04S4G", 
  "E06S1D", "E06S4D", "E06S3D", 
  "E06S1G", "E06S4G", "E06S3G"
)

gsva_res_filtered <- gsva_res_filtered[sample_order, ]

### Step 4: Prepare Heatmap Annotation ###
# Assume `column_annotation` is a properly defined HeatmapAnnotation object
# If not defined yet, you should create it like this (example):
# column_annotation <- HeatmapAnnotation(
#   Treatment = sample_conditions_vector,
#   col = list(Treatment = condition_colors)
# )

### Step 5: Plot GSVA Heatmap ###
ht <- Heatmap(
  t(gsva_res_filtered),  # Pathways as rows, samples as columns
  top_annotation = column_annotation,
  name = "GSVA score",
  cluster_columns = FALSE,
  heatmap_legend_param = list(
    direction = "vertical"
  )
)

# Draw the heatmap with adjusted padding
draw(
  ht, 
  heatmap_legend_side = "left", 
  annotation_legend_side = "left",
  padding = unit(c(1, 1, 1, 5), "cm")  # top, right, bottom, left
)
