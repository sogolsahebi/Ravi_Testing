# Clinical Data Processing for Ravi_dataset.

# Libraries and Source Files
# Access necessary functions from ICB_common codes.
source("https://raw.githubusercontent.com/BHKLAB-Pachyderm/ICB_Common/main/code/Get_Response.R")
source("https://raw.githubusercontent.com/BHKLAB-Pachyderm/ICB_Common/main/code/format_clin_data.R")
source("https://raw.githubusercontent.com/BHKLAB-Pachyderm/ICB_Common/main/code/annotate_tissue.R")


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

clin_merge_data <- read.table( file_path , stringsAsFactors=FALSE , sep="\t" )

# Data Curating
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
  "WES_All",
  "RNA_All"
)


# Combine selected columns with additional columns.
clin <- cbind(clin_merge_data[, selected_cols], "Lung", NA , NA, NA, NA, NA, NA)

# Set new column names.
colnames(clin) <- c(
  "patient", "age", "recist", "t.os", "t.pfs", "histo", "stage", "sex", "os",
  "pfs", "dna", "rna", "primary", "drug_type", "response.other.info", "response"
)


# Modify stage values for better clarity.
clin$stage <- ifelse(clin$stage == 1, "I",
                     ifelse(clin$stage == 2, "II",
                            ifelse(clin$stage == 3, "III",
                                   ifelse(clin$stage == 4, "IV", NA))))

# Assign 'tpm' and 'wes' values based on RNA_All and WES_All respectively.
clin$rna <- ifelse(clin_merge_data$RNA_All == 1, "tpm", NA)
clin$dna <- NA

# Calculate the response using Get_Response function.
clin$response = Get_Response(data = clin)


# Reorder columns.
clin <- clin[, c(
  "patient", "sex", "age", "primary", "histo", "stage", 
  "response.other.info", "recist", "response", "drug_type", "dna", "rna", "t.pfs", 
  "pfs", "t.os", "os"
)]


# Use the format_clin_data function for further formatting.
clin <- format_clin_data(clin_merge_data, "patient", selected_cols, clin)

#add survival_unit and survival_type columns. and use  'curation_tissue.csv' file to set clin$tissueid column
#annotation_tissue <- read.csv("~/BHK lab/Ravi/deleted repo/Ravi_version2/Common files/curation_tissue.csv")

#Read curation_tissue.cv file

path <- "https://raw.githubusercontent.com/BHKLAB-DataProcessing/ICB_Common/main/data/curation_tissue.csv"
#annotation_tissue<- read.csv(file.path(annot_dir, 'curation_tissue.csv'))
annotation_tissue<- read.csv(path )

clin <- annotate_tissue(clin=clin, study='Ravi', annotation_tissue=annotation_tissue, check_histo=FALSE)

#add column treatment id after tissueid column.
clin <- add_column(clin, treatmentid= clin$Agent_PD1, .after='tissueid')

# Adding drug_type based on treatmentid.
# Print unique values of treatmentid.
print(unique(clin$treatmentid))

clin$drug_type[clin$treatmentid %in% c('Nivolumab', 'Pembrolizumab', 'Atezolizumab', 'Avelumab')] <- 'PD-1/PD-L1'
clin$drug_type[clin$treatmentid %in% c('Nivolumab + LAG-3', 'Atezolizumab + Epacadostat', 'Nivolumab + Urelumab',
                                       'Nivolumab + Lirilumab', 'Durvalumab + Tremelimumab', 'Nivolumab + Ipilimumab')] <- 'IO+combo'
clin$drug_type[clin$treatmentid %in% c('Pembrolizumab + Carboplatin + Pemetrexed', 'Nivolumab + Gemcitabine')] <- 'IO+chemo'

# Print unique values of drugtype
print(unique(clin$drug_type))

# Replace empty string values with NA.
clin[clin == ""] <- NA

# Replace hyphens in the row names of 'clin' with underscores
######################## comment: clin$patient needs to be fixed as well, ow not matched samples' name between clin and expr data #########################
clin$patient <- str_replace_all(clin$patient, '-', '_') #? not sure

# Save the processed data as CLIN.csv file
#write.csv(clin , file=file.path(output_dir, "CLIN.csv") , row.names=TRUE )

file <- "~/BHK lab/Ravi_Testing/files/CLIN.csv"
write.csv(clin, file, row.names = TRUE)



