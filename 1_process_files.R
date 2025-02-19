# read some files
library(readxl)
abundance <- read.csv("C:/Users/Orian/OneDrive/Documents/University/Master/THESIS/data/sequencing_data/pathabundance.csv", header = TRUE, sep = "\t")
metadata <- read_excel("C:/Users/Orian/OneDrive/Documents/University/Master/THESIS/data/OSM metadata smokers & contact lenses.xlsx")               

# first, filter abundance based on metadata to remove smoker --> remove all columns that have CR in the name
# keep only abundance
abundance_lense <- abundance[, !grepl("CR", colnames(abundance))]
abundance_lense <- abundance_lense[, !grepl("H150", colnames(abundance_lense))]
abundance_lense <- abundance_lense[, !grepl("Abundance", colnames(abundance_lense))]
# also remove negative and posisitve control
abundance_lense <- abundance_lense[, !grepl("Pos", colnames(abundance_lense))]
abundance_lense <- abundance_lense[, !grepl("Neg", colnames(abundance_lense))]
abundance_lense <- abundance_lense[, !grepl("Plus", colnames(abundance_lense))]
# remove pathways since we are only interested in abundances
abundance_lense <- abundance_lense[!grepl("PWY", abundance_lense$X..Pathway), ]
abundance_lense <- abundance_lense[!grepl("HEME", abundance_lense$X..Pathway), ]
abundance_lense <- abundance_lense[!grepl("GLYCOLYSIS", abundance_lense$X..Pathway), ]



write.csv(abundance_lense, "abundance_lense.csv", row.names = F)

# modify metadata to match abundance
metadata_lense <- metadata[!grepl("CR", metadata$sample_name), ]
metadata_lense <- metadata_lense[!grepl("Pos", metadata_lense$sample_name), ]
metadata_lense <- metadata_lense[!grepl("Neg", metadata_lense$sample_name), ]
metadata_lense$sample <- gsub("_$", "", metadata_lense$sample)
metadata_lense$sample <- paste0("X", metadata_lense$sample) # need to manually remove the X before the B
write.csv(metadata_lense, "metadata_lense.csv", row.names = F)


