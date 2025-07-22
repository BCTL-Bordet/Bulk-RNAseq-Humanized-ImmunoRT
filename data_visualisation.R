##############################################
# Sample Clustering and PCA
# Description: 
#   - Perform hierarchical clustering using sample distances
#   - Visualize dendrogram with treatment group coloring
#   - Run PCA on log-transformed TPM matrix
#   - Generate PCA plot with treatment group coloring
##############################################

### Load Required Libraries ###
library(dendextend)
library(ggplot2)

### Step 1: Hierarchical Clustering ###
# Create sample group labels (3 replicates per condition)
classes <- rep(
  c("Saline", 
    "3x8Gy + Saline", 
    "3x8Gy + Saline dist", 
    "3x8Gy + Pembro", 
    "3x8Gy + Pembro dist"), 
  each = 3
)
names(classes) <- colnames(dist_matrix)

# Define custom colors for each group
class_colors <- c(
  "Saline" = "lightblue",
  "3x8Gy + Saline" = "yellow3",
  "3x8Gy + Saline dist" = "purple4",
  "3x8Gy + Pembro" = "green4",
  "3x8Gy + Pembro dist" = "orange3"
)

# Compute hierarchical clustering
hc <- hclust(dist_obj)
dend <- as.dendrogram(hc)

# Plot dendrogram
plot(
  dend, 
  main = "Hierarchical Clustering Dendrogram",
  xlab = "Samples", ylab = "Distance"
)
legend(
  "topright", 
  legend = names(class_colors),
  col = class_colors, 
  pch = 16,
  title = "Treatment group", 
  bty = "n",
  cex = 0.7,
  x.intersp = 0.6, 
  y.intersp = 0.6,
  inset = c(-0.55, -0.2)
)

### Step 2: PCA on log-transformed TPM ###
# Transpose matrix to have samples as rows
tpm_t <- t(log_transformed_tpm)

# Perform PCA with scaling
pca_res <- prcomp(tpm_t, scale. = TRUE)

# Extract percent variance explained by PC1 and PC2
percentVar <- round(100 * (pca_res$sdev^2 / sum(pca_res$sdev^2)), 1)

# Create data frame for ggplot
pca_data <- as.data.frame(pca_res$x[, 1:2])
pca_data$sample <- rownames(pca_data)
pca_data$condition <- rep(
  c("Saline", 
    "3x8Gy + Saline", 
    "3x8Gy + Saline dist", 
    "3x8Gy + Pembro", 
    "3x8Gy + Pembro dist"), 
  each = 3
)

# Ensure consistent ordering in plot
pca_data$condition <- factor(
  pca_data$condition,
  levels = c("Saline", "3x8Gy + Saline", 
             "3x8Gy + Saline dist", 
             "3x8Gy + Pembro", 
             "3x8Gy + Pembro dist")
)

# Define custom color palette
condition_colors <- c(
  "Saline" = "lightblue",
  "3x8Gy + Saline" = "yellow3",
  "3x8Gy + Saline dist" = "purple4",
  "3x8Gy + Pembro" = "green4",
  "3x8Gy + Pembro dist" = "orange3"
)

# Generate PCA plot
ggplot(pca_data, aes(x = PC1, y = PC2, color = condition, label = sample)) +
  geom_point(size = 4, alpha = 0.8) +
  scale_color_manual(values = condition_colors) +
  theme_minimal() +
  labs(
    title = "PCA Plot (log2 TPM)", 
    x = paste0("PC1 (", percentVar[1], "%)"), 
    y = paste0("PC2 (", percentVar[2], "%)"),
    color = "Treatment group"
  ) +
  theme(
    legend.title = element_text(face = "bold")
  )
