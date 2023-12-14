# RNA-seq Data Processing.
# File: Format_EXPR.R

#Read libararies.
library(data.table)

# Parse command line arguments for input, output, and annotation directory paths.
#args <- commandArgs(trailingOnly = TRUE)
#input_dir <- args[1]
#output_dir <- args[2]

# Data Reading
# Define the path to open EXPR.txt.gz file.
# Read the RNA-Seq data from the gct file.
#expr <- read.csv(file.path(input_dir, "EXPR.txt.gz"), stringsAsFactors=FALSE , sep="\t" )
path <- "~/BHK lab/Ravi_Testing/files/EXPR.txt.gz"
expr <- read.csv(path, stringsAsFactors = FALSE, sep = "\t", check.names = FALSE)
colnames(expr) 

# Data Filtering
# Define the path for the 'cased_sequenced.csv' file
#case = read.csv( file.path(output_dir, "cased_sequenced.csv") , sep=";" )

file_path <- "~/BHK lab/Ravi_Testing/files/cased_sequenced.csv"
# Read the 'case' dataset
case <- read.csv(file_path, sep = ";")

# Filter the 'expr' dataset to include only patients with expr value of 1 in the 'case' dataset
expr <- expr[ , case[case$expr %in% 1,]$patient]

# Check the range of data values
range(expr) #0 214022

# Data Transformation
# Convert TPM data to log2-TPM for consistency with other data formats
expr <- log2(expr + 0.001)

dim(expr)
# Check the updated range of data values
range(expr) #now  -9.965784 17.707400

# Data Export
# Define the output path for the cleaned data.
file_path <- "~/BHK lab/Ravi_Testing/files/EXPR.csv"
write.csv(expr, file_path, row.names = TRUE)

# Write the  'expr' dataset to a CSV file.
#write.table( expr , file=file.path(output_dir, "EXPR.csv") , row.names=TRUE )



