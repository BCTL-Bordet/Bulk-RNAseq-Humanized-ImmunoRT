#### Gene Set Variation Analysis (GSVA) ####

library(GSVA)
library(ComplexHeatmap)

# Load gene sets from GMT file
gene_sets <- gmtPathways("your/pathways/file.gmt")

# Remove unwanted pathways
gene_sets_filtered <- gene_sets[!names(gene_sets) %in% c(
  "HALLMARK_PANCREAS_BETA_CELLS", "HALLMARK_SPERMATOGENESIS", 
  "HALLMARK_UV_RESPONSE_DN", "HALLMARK_UV_RESPONSE_UP", 
  "HALLMARK_ANDROGEN_RESPONSE", "HALLMARK_APICAL_SURFACE"
)]

# Perform GSVA analysis
params <- gsvaParam(as.matrix(log_transformed_tpm),   # Gene expression matrix
                    kcdf="Gaussian", geneSets = gene_sets_filtered)

# Convert results to a dataframe
gsva_res_filtered <- as.data.frame(t(gsva(params)))

# Format pathway names for readability
gsub_names <- function(names) {
  names <- gsub("HALLMARK_", "", names)
  names <- gsub("_", " ", names)
  return(names)
}
colnames(gsva_res_filtered) <- gsub_names(colnames(gsva_res_filtered))

# Define sample names for ordering
sample_names <- c("E01S1G", "E01S3G", "E01S5G", "E02S2G", "E02S3G", "E02S4G", 
                  "E03S1D", "E03S4D", "E04S4D", "E03S1G", "E03S4G", "E04S4G", 
                  "E06S1D", "E06S4D", "E06S3D", "E06S1G", "E06S4G", "E06S3G")

# Reorder results based on sample names
gsva_res_filtered <- gsva_res_filtered[sample_names, ]

# Create heatmap annotation (assumed `column_annotation` is predefined)
ht <- Heatmap(
  t(gsva_res_filtered),
  top_annotation = column_annotation,
  name = "GSVA score",
  cluster_columns = FALSE,
  heatmap_legend_param = list(direction = "vertical")
)

# Draw heatmap with custom margins
draw(
  ht, 
  heatmap_legend_side = "left", 
  annotation_legend_side = "left",
  padding = unit(c(1, 1, 1, 5), "cm")  # Top, right, bottom, left margins
)

