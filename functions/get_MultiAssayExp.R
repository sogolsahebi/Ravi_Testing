library(Biobase)
library(SummarizedExperiment)
library(MultiAssayExperiment)
library(stringr)
#add
library(readr)


get_MultiAssayExp <- function(study, expr_with_counts_isoforms=FALSE){ 

  
  library(readr)
  path <- "https://raw.githubusercontent.com/BHKLAB-DataProcessing/ICB_Common/main/data/DATASET_LOAD_INFO.csv"
  DATASET_LOAD_INFO <- read_delim("files/cased_sequenced.csv", delim = ";", escape_double = FALSE, trim_ws = TRUE)
  
  
  # Create a new row for Ravi 
  #new_row <- data.frame(
  #  study = "Ravi", 
  #  expr_bool = TRUE, 
  #  snv_bool = FALSE,  
  #  cna_bool = FALSE, 
  #  cin_bool = FALSE, 
  # coverage = NA, 
  #  indel_bool = FALSE, 
  #  mutsig_bool = FALSE, 
  # WES.comment = NA,
  # stringsAsFactors = FALSE)
  
  #for Ravi we have :
  study = "Ravi" 

  
  #add input_diras alist 
  input_dir2 = list()
  input_dir2[["CLIN"]] <- read.csv("files/CLIN.csv",row.names = 1)
  input_dir2[["EXPR"]] <- read.csv("files/EXPR.csv",row.names = 1)
  input_dir2[["case"]] <- read_delim("files/cased_sequenced.csv", delim = ";", escape_double = FALSE, trim_ws = TRUE)
  
  
  #se_list <- Create_SummarizedExperiments( 
    #input_dir= input_dir2,
    #study= data$study, 
    #expr_bool= data$expr_bool, 
    #snv_bool= data$snv_bool, 
    #cna_bool= data$cna_bool, 
    #cin_bool= data$cin_bool, 
    #coverage= data$coverage, 
    #indel_bool= data$indel_bool,
    #expr_with_counts_isoforms=expr_with_counts_isoforms
  #)
  
  
  #add
  se_list <- list()
  
  se_list[["expr"]] <- Create_EXP_SummarizedExperiment(
    study=Ravi , 
    case=case, 
    clin=clin, 
    expr=expr, 
    feat_snv=FALSE , 
    feat_cna=FALSE, 
    feat_cin=FALSE, 
    cna_bool=FALSE , 
    snv_bool=FALSE 
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
  #ICB_Ravi <- MultiAssayExperiment(experiments=se_list, colData=coldata)

  # Save the multiassay_result object as an RDS file
  #saveRDS(ICB_Ravi, file = "C:/Users/sogol/OneDrive/Documents/BHK lab/Ravi_version2/data/ICB_Ravii.rds")

  

  return(MultiAssayExperiment(experiments=se_list, colData=coldata))
}

