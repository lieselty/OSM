# Adapted from Nina Eldridge

# Load library
library(ggplot2)
library(vegan)
library(dplyr)
library(tidyr)

# read in the data and the metadata
data <- read.csv("C:/Users/Orian/OneDrive/Documents/University/Master/THESIS/scripts/abundance_lense.csv", header = TRUE)
rownames(data) <- data$X..Pathway
data <- data[,-1]


# set a seed for result reproducability
set.seed(123)


###### basic statistics #######
# orientation: columns = samples, rows = species

# gets you the mean of each row
row_means <- rowMeans(data, na.rm = TRUE)
# gets you the median of each row
row_medians <- apply(data, 1, median, na.rm = TRUE)

# gets you the proportion of how many samples the species is present
presence_proportion <- apply(data, 1, function(row) {
  # Count non-zero values and divide by the number of columns
  sum(row != 0, na.rm = TRUE) / (ncol(data))
})

# gets you the amount of samples where the species was detected
sample_counts <- apply(data, 1, function(row) {
  # Count non-zero values and divide by the number of columns
  sum(row != 0, na.rm = TRUE)
})

data_info <- data
data_info$means <- row_means
data_info$medians <- row_medians
data_info$sample_counts <- sample_counts
data_info$presence_proportion <- presence_proportion

###### alpha diversity ######
shannon_index <- diversity(data, index = "shannon", MARGIN = 2) # high value = high diversity
simpson_index <- diversity(data, index = "simpson", MARGIN = 2) # high values = low density
head(shannon_index)
head(simpson_index)

###### plots #######

### cumulative plot ###
# calculate the prevalence at each threshold
count_prevalence_steps <- function(column) {
  # Define the sequence of prevalence thresholds
  # Adapt the value until where, such that we can see small counts if necessary
  thresholds <- seq(1, 0, by = -0.05)
  
  # Count the occurrences for each threshold
  counts <- sapply(thresholds, function(threshold) sum(column >= threshold))
  
  # Return the counts as a named vector
  names(counts) <- thresholds
  return(counts)
}

# plots the calculated prevalences for each threshold
plot_prevalence_counts <- function(prevalence_counts) {
  # Plot the results
  plot(
    x = names(prevalence_counts), # Names of the vector as x-axis labels
    y = prevalence_counts,        # Counts as y values
    type = "o",                   # "o" for lines and points
    xlab = "Prevalence Threshold", 
    ylab = "Counts", 
    main = "Number of species with a certain prevalence",
    col = "violet",
    pch = 16 # Filled circles for points
  )
  
  # Rotate x-axis labels for better readability
  axis(1, at = seq_along(names(prevalence_counts)), labels = names(prevalence_counts), las = 2, cex.axis=0.8)
  
  # Add grid lines for better readability
  grid()
}

# example usage
result <- count_prevalence_steps(data_info$presence_proportion)
plot_prevalence_counts(result)


### composition ###
# Filter the data to remove rows and columns with all zeros
data_compo <- data[rowSums(data != 0) > 0, ]
data_compo <- data_compo[, colSums(data_compo != 0) > 0]

# change the values of rows with UNINTEGRATED to have the correct fraction (I think before it was a fraction of UNINTEGRATED but I want the total fraction)
UNINTEGRATED_row <- as.matrix(data_compo["UNINTEGRATED", , drop = FALSE])
rows_to_update <- grepl("^UNINTEGRATED", rownames(data_compo)) & rownames(data_compo) != "UNINTEGRATED"
data_compo[rows_to_update, ] <- as.matrix(data_compo[rows_to_update, ]) * UNINTEGRATED_row[rep(1, sum(rows_to_update)), ]
sum_UNINTEGRATED_rows <- colSums(data_compo[rows_to_update, , drop = FALSE])
data_compo["UNINTEGRATED", ] <- UNINTEGRATED_row - sum_UNINTEGRATED_rows
# add an unkown column
unknown_values <- 1 - colSums(data_compo)
data_compo <- rbind(data_compo, "0_UNKNOWN" = unknown_values)

# Plot the composition
data_compo$Taxonomic_Name <- rownames(data_compo)

df_long <- data_compo %>%
  pivot_longer(cols = -Taxonomic_Name, names_to = "Sample", values_to = "Abundance")

color_palette <- c("0_UNKNOWN" = "grey80", "UNMAPPED" = "gray40", "UNINTEGRATED" = "lightsalmon", 
                   "UNINTEGRATED|g__Comamonas.s__Comamonas_aquatica" = "skyblue1", 
                   "UNINTEGRATED|g__Enhydrobacter.s__Enhydrobacter_aerosaccus" = "palegreen2",
                   "UNINTEGRATED|g__Moraxella.s__Moraxella_osloensis" = "deeppink2",
                   "UNINTEGRATED|g__Staphylococcus.s__Staphylococcus_aureus" = "lightgoldenrod1",
                   "UNINTEGRATED|unclassified" = "grey60") # Add other taxa colors as needed

ggplot(df_long, aes(x = Sample, y = Abundance, fill = Taxonomic_Name)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(x = "Sample", y = "Abundance", fill = "Taxonomic Name") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_manual(values = color_palette)


### plot alpha diversity ###
# Bar plot for Shannon Index with customized x-axis labels
barplot(shannon_index, 
        main = "Shannon Diversity Index", 
        xlab = "Samples", 
        ylab = "Shannon Index", 
        col = "violet",
        las = 2,          
        cex.axis = 0.8,   
        cex.names = 0.6) 

# Bar plot for Simpson Index with customized x-axis labels
barplot(simpson_index, 
        main = "Simpson Diversity Index", 
        xlab = "Samples", 
        ylab = "Simpson Index", 
        col = "purple",
        las = 2,          
        cex.axis = 0.8,   
        cex.names = 0.6)

