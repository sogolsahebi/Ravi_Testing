library(Biobase)
library(SummarizedExperiment)
library(MultiAssayExperiment)
library(stringr)

source_location <- "https://raw.githubusercontent.com/BHKLAB-Pachyderm/ICB_Common/main"
# source_location <- "~/Documents/GitHub/Pachyderm/PredictIO_ICB/ICB_Common"

get_MultiAssayExp <- function(study, input_dir, expr_with_counts_isoforms=FALSE){ 

  #add
  library(readr)
  DATASET_LOAD_INFO <- read_csv("data/DATASET_LOAD_INFO.csv")
  
  new_row <- data.frame(
    study = "Ravi", 
    expr_bool = TRUE, 
    snv_bool = FALSE,  
    cna_bool = FALSE, 
    cin_bool = FALSE, 
    coverage = NA, 
    indel_bool = FALSE, 
    mutsig_bool = FALSE, 
    WES.comment = NA,
    stringsAsFactors = FALSE
  )
  
  
  DATASET_LOAD_INFO <- rbind(DATASET_LOAD_INFO, new_row)
  data <- data[data$study == study, ]
  
  
  
  se_list <- Create_SummarizedExperiments( 
    input_dir=input_dir,
    study= data$study, 
    expr_bool= data$expr_bool, 
    snv_bool= data$snv_bool, 
    cna_bool= data$cna_bool, 
    cin_bool= data$cin_bool, 
    coverage= data$coverage, 
    indel_bool= data$indel_bool,
    expr_with_counts_isoforms=expr_with_counts_isoforms
  )
  
  
  cols <-list()
  for(assay_name in names(se_list)){
    cols[[assay_name]]<- data.frame(colData(se_list[[assay_name]]))
  }
  
  #cols[["expr"]] <- coldata #152 
  
  
  # Format and merge coldata
  allcols <- lapply(cols, function(col){
    return(rownames(col))
  })


  allcols <- unique(unlist(allcols))
  #length(allcols) 152
  
  coldata <- NA
  
  for(col in cols){
    if(!is.data.frame(coldata)){
      coldata <- col
    }
    missing <- allcols[!allcols %in% rownames(coldata)]
    filtered <- col[rownames(col) %in% missing, ]
    if(nrow(filtered) > 0){
      coldata <- rbind(coldata, filtered)
    }
  }
  
  coldata <- coldata[order(rownames(coldata)), ]
  dim(coldata)
  
  #add
  ICB_Ravi <- MultiAssayExperiment(experiments=se_list, colData=coldata)

  # Save the multiassay_result object as an RDS file
  saveRDS(ICB_Ravi, file = "C:/Users/sogol/OneDrive/Documents/BHK lab/Ravi_version2/data/ICB_Ravii.rds")

  e <- assays(ICB_Ravi)[["expr"]]
  c <- data.frame(colData(ICB_Ravi))
  a <- data.frame(rowData(ICB_Ravi@ExperimentList$expr)) # 2794 protein-coding genes???? 
  
  dim(e)
  dim(c)
  dim(a)
  View(a)
  protein_coding_gene <- a[a$gene_type == "protein_coding",]
  dim(protein_coding_gene)
  
  return(MultiAssayExperiment(experiments=se_list, colData=coldata))
}
