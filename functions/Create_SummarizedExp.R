library(data.table)
library(Biobase)
library(SummarizedExperiment)
library(GenomicRanges)
library(biomaRt)
library(stringr)


clin_cols <- c(
  "patient" , "sex" , "age" , "primary" , "histo" , "tissueid", "treatmentid", "stage" , 
  "response.other.info" , "recist" , "response" , "drug_type" , 
  "dna" , "rna" , "t.pfs" , "pfs" , "t.os" , "os", 
  "survival_unit", "survival_type"
)

added_cols <- c(
  "TMB_raw" , "nsTMB_raw" , "indel_TMB_raw" , "indel_nsTMB_raw" , 
  "TMB_perMb" , "nsTMB_perMb" , "indel_TMB_perMb" , "indel_nsTMB_perMb" ,
  "CIN" , "CNA_tot" , "AMP" , "DEL" 
)

renamed_cols <- list(
  drug_type = "treatment",
  primary = "cancer_type",
  t.pfs = "survival_time_pfs",
  pfs = "event_occurred_pfs",
  t.os = "survival_time_os",
  os = "event_occurred_os",
  patient = "patientid"
)

format_se <- function(assay, coldata, assay_type, convert_gene_name=TRUE, is_isoform=FALSE){
  
  
  #add
  colnames(expr) <- str_replace_all(colnames(expr), '[-\\.]', '_')
  rownames(clin) <- str_replace_all(rownames(clin), '[-\\.]', '_')
  patient = intersect( colnames(expr) , rownames(clin) )
  colnames(expr)
  rownames(clin)
  patient
  coldata = clin[ patient , ]
  
  
  assay = expr
  assay_type = "expr"
  
  for(renamed_col in names(renamed_cols)){
    colnames(coldata)[colnames(coldata) == renamed_col] <- renamed_cols[[renamed_col]]
  }
  
  
  #loag genecode version 19 
  load("~/BHK lab/Annotation/Gencode.v19.annotation (2).RData")
  gene_ids <- c()
  features_df <- features_gene
  
  
  ##add Sisira 
  rownames(features_df) <- gsub("\\..*","", rownames(features_df))
  rownames(features_df)
  
  
  #add
  assay$ensembel_noversion <- gsub("\\..*","", rownames(assay))
  
  
  #length(intersect(rownames(features_df), assay$ensembel_noversion))
  #add(changed)
  assay <- assay[assay$ensembel_noversion %in% rownames(features_df), ]
  assay <- assay[!duplicated(assay$ensembel_noversion), ]
  dim(assay)


  #add
  #assay$ensembel_noversion <- NULL
  
  gene_ids <- rownames(assay)
  

  #add assay$ensembel_noversion
  #build the GRanges object ussed as rowRanges (rowData)
  #assay_genes <- as.data.table(features_df[rownames(features_df) %in% assay$ensembel_noversion, ])
  #change to data.table 
  assay_genes <- data.frame(features_df[rownames(features_df) %in% assay$ensembel_noversion, ])

  
  #chcek duplicate
  duplicates_in_features_df <- any(duplicated(rownames(features_df)))
  duplicates_in_assay <- any(duplicated(assay$ensembel_noversion))
  # Count the number of duplicated row names in features_df
  num_duplicates_in_assay <- sum(duplicated(assay$ensembel_noversion))
  num_duplicates_in_assay
  duplicates_in_features_df
  duplicates_in_assay
  
  
  dim(assay_genes)
  rownames(assay_genes)
  
  #add (maybe redunanatnt since later is coneverting to  data table) data,table??
  rownames(assay_genes) <- assay_genes$gene_id
  rownames(assay_genes)
  
  library(data.table)
  assay_genes <- assay_genes[order(rownames(assay_genes)), ]
  
  #check this
  #assay_genes[, gene_id_no_ver := gsub("\\..*$", "", gene_id)]
  
  #assay_genes[
  #  is.na(start),
  #  c("start", "end", "length", "strand") := list(-1, -1, 0, "*")
  #]
  

  assay_genes <- assay_genes[order(assay_genes$gene_id), ]
  
  row_ranges <- makeGRangesFromDataFrame(
    assay_genes,
    keep.extra.columns=TRUE  # retain metadata
  )
  
  names(row_ranges) <- row_ranges$rownames
  
  assay_list <- list()
  #add
  assay$ensembel_noversion <- NULL
  
  assay_list[[assay_type]] <- assay


  #add
  dim(assay)
  summary(row_ranges)
  dim(assay_list[["expr"]])
  dim(coldata)

  length(row_ranges)
  
  S1 <- SummarizedExperiment(assays=assay_list, colData=coldata, rowRanges=row_ranges)
  
  
  
  return(SummarizedExperiment(assays=assay_list, colData=coldata, rowRanges=row_ranges))

  
}



Create_EXP_SummarizedExperiment = function( study, case , clin, expr, feat_snv, feat_cna, feat_cin, snv_bool, cna_bool, is_isoform=FALSE ){
  
  

  rownames(clin) = clin$patient
  added_df <- as.data.frame(matrix(NA, nrow = nrow(clin), ncol = length(added_cols)))
  colnames(added_df) <- added_cols
  clin <- data.frame(cbind(
    clin[, clin_cols],
    added_df,
    clin[, !colnames(clin) %in% clin_cols]
  ))
  
  #add
  rownames(clin)
  rownames(clin) <- str_replace_all(rownames(clin), '[-\\.]', '_')
  colnames(expr) <- str_replace_all(colnames(expr), '[-\\.]', '_')
  patient = intersect( colnames(expr) , rownames(clin) )
  patient
  colnames(expr)
  clin = clin[ patient , ]
  expr = expr[ , patient ]

  
  #add add is_isoform to FALSE and conver_gene_Name = FLASE
  
  return(format_se(assay=expr, coldata = clin , assay_type='expr', convert_gene_name= FALSE, is_isoform=FALSE))
}


Create_SummarizedExperiments = function( study , input_dir, expr_bool , snv_bool , cna_bool, cin_bool , coverage , indel_bool, expr_with_counts_isoforms ){
  
  
  # Path to processed data 
  #case_file = paste( input_dir , "cased_sequenced.csv" , sep="/" )
  #clin_file = paste( input_dir , "CLIN.csv" , sep="/" )
  #expr_file = paste( input_dir , "EXPR.csv" , sep="/" )
  #snv_file = paste( input_dir , "SNV.csv" , sep="/" )
  #cna_file = paste( input_dir , "CNA_gene.csv" , sep="/" )
  
  feat_snv <- NA
  feat_cna <- NA
  feat_cin <- NA
  
  se_list <- list()

  case = read_delim("files/cased_sequenced.csv", delim = ";", escape_double = FALSE, trim_ws = TRUE)
  clin = read.csv("files/CLIN.csv",row.names = 1,check.names = FALSE)
  expr = read.csv("files/EXPR.csv",row.names = 1, check.names = FALSE)
  
  if( expr_bool ){
    
    
    colnames(expr) <- str_replace_all(colnames(expr), '[-\\.]', '_')
    
    se_list[["expr"]] <- Create_EXP_SummarizedExperiment(
      study=study, 
      case=case, 
      clin=clin, 
      expr=expr, 
      feat_snv=feat_snv , 
      feat_cna=feat_cna, 
      feat_cin=feat_cin , 
      cna_bool=cna_bool , 
      snv_bool=snv_bool 
    ) 
    
  }
  
  return(se_list)
}

