# Format_downloaded_data.R.

# This script formats and cleans clinical and expression data. 
# - Creates "CLIN.txt" by merging clinical data from two sources.
# - Creates "EXPR.txt.gz".


# Load necessary libraries
library(readxl)
library(CePa)

# Setup working directory from command line arguments
#args <- commandArgs(trailingOnly = TRUE)
#work_dir <- args[1]

# Alternatively, set file paths directly
#file_path1 <- file.path(work_dir,"SU2C-MARK_Harmonized_Clinical_Annotations_Supplement_v1.txt" )
#file_path2 <- file.path(work_dir,"Table_S1_Clinical_Annotations.xlsx" )


file_path1 <- "~/BHK lab/Ravi_Testing/Source Data/Source Data/Clinical/SU2C-MARK_Harmonized_Clinical_Annotations_Supplement_v1.txt"
file_path2 <- "~/BHK lab/Ravi/Source Data/Source Data/Clinical/Table_S1_Clinical_Annotations.xlsx"

# STEP 1: Create "CLIN.txt" file
# Load clinical data

clin_data1 <- read.table(file_path1, header = TRUE, sep = "\t")
clin_data2 <- read_excel(file_path2, col_names = FALSE)


# Clean clinical data, Set row 3 as column names  
colnames(clin_data2) <- clin_data2[3, ]
# Remove the first 3 rows
clin_data2 <- clin_data2[-(1:3), ]


# Merge clinical data
clin_merge_data <- merge(clin_data1, clin_data2, by = "Harmonized_SU2C_Participant_ID_v2")


# Change the first column name to patient
colnames(clin_merge_data)[colnames(clin_merge_data) == "Harmonized_SU2C_Participant_ID_v2"] <- "patient"



# Save merged clinical data
# write.table(clin_merge_data, file=file.path(work_dir, 'CLIN.txt'), quote=FALSE , sep="\t" , col.names=TRUE , row.names=FALSE)

file_path <- "~/BHK lab/Ravi_Testing/files/CLIN.txt"
write.table(clin_merge_data, file_path , sep = "\t", row.names = TRUE)


# STEP 2: Create "EXPR.txt.gz" file
# Define the path for the .gct file.
#gct_file_path <- file.path(work_dir,"SU2C-MARK_Harmonized_rnaseqc_tpm_v1.gct" )

gct_file_path<- "~/BHK lab/Ravi_Testing/Source Data/Source Data/RNA/SU2C-MARK_Harmonized_rnaseqc_tpm_v1.gct"

# Load expression data
expr <- data.frame(read.gct(gct_file_path))

# Clean expression data
# 1. Replace periods with hyphens
# 2. Remove trailing -T1 or -T2 from colnames(expr)
new_colnames <- gsub("\\.", "-", colnames(expr))
new_colnames <- gsub("-T1$|-T2$", "", new_colnames)
colnames(expr) <- new_colnames


# Check for any duplicate column names after renaming
duplicate_colnames <- colnames(expr)[duplicated(colnames(expr))]
print(duplicate_colnames)  #should print nothing ideally

# Sort the row names of 'expr'
expr <- expr[sort(rownames(expr)),]


# Confirm that all column values are numeric
stopifnot(all(sapply(expr, is.numeric)))


# Open a gzipped file for writing
#gz <- gzfile(file.path(work_dir, 'EXPR.txt.gz'), "w")

gz <- gzfile("~/BHK lab/Ravi_Testing/files/EXPR.txt.gz", "w")

# Use write.table to write the data to the gzipped file
write.table(expr, file = gz, sep = "\t", row.names = TRUE, quote = FALSE)

# Close the gzipped file
close(gz)

