# Clinical Data Processing
# Goal: save CLIN.csv (dimensions: 389 x 54).

# Libraries and Source Files
# Access necessary functions from ICB_common codes.
source("https://raw.githubusercontent.com/BHKLAB-Pachyderm/ICB_Common/main/code/Get_Response.R")
source("https://raw.githubusercontent.com/BHKLAB-Pachyderm/ICB_Common/main/code/format_clin_data.R")
source("https://raw.githubusercontent.com/BHKLAB-Pachyderm/ICB_Common/main/code/annotate_tissue.R")
source("https://raw.githubusercontent.com/BHKLAB-Pachyderm/ICB_Common/main/code/annotate_drug.R")

#read library tibble for add_column function
library(tibble)

# Parse command line arguments for input, output, and annotation directory paths.
#args <- commandArgs(trailingOnly = TRUE)
#input_dir <- args[1]
#output_dir <- args[2]
#annot_dir <- args[3]

# Data Loading
#Define the path to the source data.
file_path <- "~/BHK lab/Ravi_Testing/files/CLIN.txt"

# Load the clinical merged data from the given file path.
#clin_merge_data <- read.csv( file.path(input_dir, "CLIN.txt") , stringsAsFactors=FALSE , sep="\t" )
clin_merged <- read.table( file_path , stringsAsFactors=FALSE , sep="\t")

# Select the required columns for further analysis.
selected_cols <- c(
  "patient",
  "Patient_Age_at_Diagnosis",
  "Harmonized_Confirmed_BOR",
  "Harmonized_OS_Months",
  "Harmonized_PFS_Months",
  "Histology_Harmonized",
  "Initial_Stage",
  "Patient_Sex",
  "Harmonized_OS_Event",
  "Harmonized_PFS_Event", 
  "Agent_PD1")

# Combine selected columns with additional columns.
clin <- cbind(clin_merged[, selected_cols], "Lung", NA , NA, NA, NA, NA ,NA )

# Set new column names
colnames(clin) <- c(
  "patient", "age", "recist", "t.os", "t.pfs", "histo", "stage", "sex", "os",
  "pfs", "drug_type", "primary", "dna", "dna_info", "rna", "rna_info", "response.other.info", "response")

# Modify stage values
clin$stage <- ifelse(clin$stage == 1, "I",
                     ifelse(clin$stage == 2, "II",
                            ifelse(clin$stage == 3, "III",
                                  ifelse(clin$stage == 4, "IV", NA))))

# Initialize DNA-related columns with NA
clin$dna <- NA
clin$dna_info <- NA

# Set RNA-related columns based on 'RNA_All' values
clin$rna <- ifelse(clin_merged$RNA_All == 1, "rnaseq", NA)
clin$rna_info <- ifelse(clin_merged$RNA_All == 1, "tpm", NA)

# Set Na for  WES related columns.
clin_merged$WES_All <- NA
clin_merged$WES_Cohort_1 <- NA
clin_merged$WES_Cohort_2 <- NA
  
# Calculate the response using Get_Response function.
clin$response = Get_Response(data = clin)

# Reorder columns.
clin <- clin[, c(
  "patient", "sex", "age", "primary", "histo", "stage", 
  "response.other.info", "recist", "response", "drug_type", "dna", "dna_info", "rna","rna_info", "t.pfs", 
  "pfs", "t.os", "os"
)]

# Use the format_clin_data function for further formatting.
clin <- format_clin_data(clin_merged, "patient", selected_cols, clin)

# Read 'curation_tissue.csv' file
# annotation_tissue<- read.csv(file.path(annot_dir, 'curation_tissue.csv'))
path <- "https://raw.githubusercontent.com/BHKLAB-DataProcessing/ICB_Common/main/data/curation_tissue.csv"
annotation_tissue<- read.csv(path)

# Annotate 'clin' using the 'annotation_tissue' data ; Set tissueid column (survival_unit and survival_type columns are added in this step).
clin <- annotate_tissue(clin=clin, study='Ravi', annotation_tissue=annotation_tissue, check_histo=FALSE) 

# Set treatmentid after tissueid column ,based on curation_drug.csv file
annotation_drug <- read.csv("https://raw.githubusercontent.com/BHKLAB-DataProcessing/ICB_Common/main/data/curation_drug.csv")
clin <- add_column(clin, treatmentid=annotate_drug('Ravi', clin$drug_type, annotation_drug), .after='tissueid')

# Check unique values of treatmentid
unique(clin$treatmentid)

# Set drug_type based on treatmentid
clin$drug_type[clin$treatmentid %in% c('Nivolumab', 'Pembrolizumab', 'Atezolizumab', 'Avelumab')] <- 'PD-1/PD-L1'
clin$drug_type[clin$treatmentid %in% c('Nivolumab + LAG-3', 'Atezolizumab + Epacadostat', 'Nivolumab + Urelumab',
                                       'Nivolumab + Lirilumab', 'Durvalumab + Tremelimumab', 'Nivolumab + Ipilimumab')] <- 'IO+combo'

clin$drug_type[clin$treatmentid %in% c('Pembrolizumab + Carboplatin + Pemetrexed', 'Nivolumab + Gemcitabine')] <- 'IO+chemo'

# Check unique values of drugtype
unique(clin$drug_type)

# Replace empty string values with NA
clin[clin == ""] <- NA

# Save the processed data as CLIN.csv file
# write.table( clin , file = file.path(output_dir, "CLIN.csv") , quote=FALSE , sep=";" , col.names=TRUE , row.names=FALSE )
file_path <- "~/BHK lab/Ravi_Testing/files/CLIN.csv"
write.table( clin , file=file_path , quote=FALSE , sep=";" , col.names=TRUE , row.names=FALSE )
