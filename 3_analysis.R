# Adapted from Nina Eldridge

# library for plotting
library(ggplot2)
# library for most analysis functions such as PCA etc
library(vegan)
# library for maaslin
library(Maaslin2)


# read in the data and the metadata
data <- read.csv("C:/Users/Orian/OneDrive/Documents/University/Master/THESIS/scripts/abundance_lense.csv", header = TRUE)
rownames(data) <- data$X..Pathway
data <- data[,-1]

metadata <- read.csv("C:/Users/Orian/OneDrive/Documents/University/Master/THESIS/scripts/metadata_lense.csv", header = T, sep= ",")               


# set a seed for result reproducability
set.seed(123)


###### change orientation of df ######
data_transposed <- data.frame(t(data))


###### analysis #####

### hclust ###
# orientation: rows = samples, columns = species

# Compute the distance matrix
dist_matrix <- dist(data_transposed)
# Perform hierarchical clustering
hc <- hclust(dist_matrix)
# Plot the dendrogram
plot(hc, main = "Dendrogram of Hierarchical Clustering", xlab = "", sub = "", cex = 0.8, labels = F)

# add colour, violent for lense 
labels <- hc$labels
label_colors <- ifelse(grepl("CK", labels) | labels == "B2OUFSH", "violet", "black")
text(x = 1:length(labels), 
     y = rep(min(hc$height) - 0.1 * diff(range(hc$height)), length(labels)), 
     labels = labels, col = label_colors, srt = 90, adj = 1, xpd = TRUE, cex = 0.8)
legend("topright", legend = c("Contact lense", "contol"), 
       col = c("violet", "black"), pch = 15, cex = 0.8, bty = "n")


### PCA ###
# orientation: rows = samples, columns = species
# Perform PCA (with scaling)
pca.result <- prcomp(data_transposed)
pca.scores <- as.data.frame(pca.result$x)

# Calculate variance explained
variance_explained <- pca.result$sdev^2
total_variance <- sum(variance_explained)
proportion_explained <- variance_explained / total_variance
percentage_explained <- proportion_explained * 100

# Print variance explained
cat("Variance explained by each PC:\n")
print(percentage_explained)

# Ensure metadata is aligned with PCA scores
metadata <- metadata[order(metadata$sample), ]
pca.scores <- pca.scores[order(rownames(pca.scores)), ]

# Bind metadata to PCA scores
pca.data <- cbind(pca.scores, treatment = metadata$treatment)

# Plot PCA results
ggplot(pca.data, aes(x = PC1, y = PC2, color = treatment)) +
  geom_point(size = 3) +
  theme_gray() +
  ggtitle("PCA contact lense") +
  xlab(paste0("Principal Component 1 [", round(percentage_explained[1], 1), "%]")) +
  ylab(paste0("Principal Component 2 [", round(percentage_explained[2], 1), "%]"))



### PCoA ###
# orientation: rows = samples, columns = species
# Step 1: Calculate a distance matrix
distance_matrix <- vegdist(data_transposed, method = "euclidian")  # Use Bray-Curtis distance --> not working, use eucilidian

# Step 2: Perform PCoA
pcoa_result <- cmdscale(distance_matrix, k = 2, eig = TRUE)

# Step 3: Create a dataframe for PCoA results
pcoa_df <- data.frame(
  PCoA1 = pcoa_result$points[, 1],
  PCoA2 = pcoa_result$points[, 2]
)

# Step 4: Combine PCoA results with metadata
pcoa_combined <- cbind(pcoa_df, metadata$treatment)
colnames(pcoa_combined)[3] <- "treatment"  # Rename the new column

# Calculate percentage of variance explained for PCoA
variance_explained <- pcoa_result$eig / sum(pcoa_result$eig) * 100

# Add percentages to axis labels
x_label <- paste0("PCoA 1 (", round(variance_explained[1], 2), "%)")
y_label <- paste0("PCoA 2 (", round(variance_explained[2], 2), "%)")

# Step 5: Plot the PCoA result
ggplot(pcoa_combined, aes(x = PCoA1, y = PCoA2, color = treatment)) +
  geom_point(size = 3) +
  labs(title = "PCoA contact lense", x = x_label, y = y_label) +
  theme_minimal() +
  theme(legend.title = element_text(size = 10),
        legend.text = element_text(size = 8))



### NMDS ###
# orientation: rows = samples, columns = species

# Step 1: Calculate a distance matrix
distance_matrix <- vegdist(data_transposed, method = "euclidian")  # Use Bray-Curtis distance

# Step 2: Perform NMDS
nmds_result <- metaMDS(distance_matrix, k = 2, trymax = 100)

# Step 3: Create a dataframe for NMDS results
nmds_df <- data.frame(
  NMDS1 = nmds_result$points[, 1],
  NMDS2 = nmds_result$points[, 2]
)

# Step 4: Combine NMDS results with metadata column you want to colour
nmds_combined <- cbind(nmds_df, metadata$treatment)
colnames(nmds_combined)[3] <- "treatment"  # Rename the new column

# Step 5: Plot the NMDS result
ggplot(nmds_combined, aes(x = NMDS1, y = NMDS2, color = treatment)) +
  geom_point(size = 3) +
  labs(title = paste0("NMDS Visualization of metaphlan OSM data"), x = "NMDS 1", y = "NMDS 2") +
  theme_minimal() +
  theme(legend.title = element_text(size = 10),
        legend.text = element_text(size = 8))


### MaAsLin ###

# run maaslin2
rownames(metadata) <- metadata$sample


fit_data <- Maaslin2(
  input_data = data,
  input_metadata = metadata,
  output = "Maaslin_output", # name of output folder
  min_abundance = 0.0,
  min_prevalence = 0.0,
  normalization = "NONE",
  transform = "LOG",
  analysis_method = "LM", # linear model
  max_significance = 0.01, # your significance level
  fixed_effects = c("treatment"), # the effect you are analysing
  random_effects = c("sex", "age"), # the difference you don't want to account for
  correction = "BH",
  standardize = FALSE,
  cores = 1)
