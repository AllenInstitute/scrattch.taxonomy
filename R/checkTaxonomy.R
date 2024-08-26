#' Checks whether an anndata object is in scrattch.taxonomy format and returns a log-file if not
#'
#' @param AIT.anndata A reference taxonomy anndata object to be tested
#' @param log.file.path The directory to output the logfile of errors and warnings (if any; default getwd())
#' 
#' Note: Any breaking issues will cause this function to return FALSE.  And potential issues will still return TRUE but will output a warning to stderr.  All messages will get returned to the log file. 
#' 
#' @import anndata
#'
#' @return Logical vector indicating whether the inputted taxonomy is a valid scrattch.taxonomy format.
#'
#' @export
checkTaxonomy = function(AIT.anndata, log.file.path=getwd()){
  
  #########################################
  ## Initial general check and set up
  
  if(sum(grepl("anndata",tolower(class(AIT.anndata))))==0){
    stop("AIT.anndata is not a variable of class anndata.")
  }
  messages  = NULL
  isValid   = TRUE
  isWarning = FALSE
  
  #########################################
  ## Check the data and metadata
  
  ## Check the normalized data (X)
  if(!is.null(AIT.anndata$X)){
    dat <- AIT.anndata$X
    if((class(dim(dat))!="integer")|(class(dat)[1]=="data.frame")){
      isValid = FALSE
      messages = c(messages,"\nERROR: LogCPM are not provided in AIT.anndata$X.")
    } else if(max(dat)>20) {
      isWarning = TRUE
      messages = c(messages,"\nWARNING: AIT.anndata$X has high values.  Please confirm this is log-normalized.")
    } else {
      messages = c(messages,":-) AIT.anndata$X looks correct.")
    }
  }
  
  ## Check the raw data (X)
  if(!is.null(AIT.anndata$raw[["X"]])){
    dat <- AIT.anndata$raw[["X"]]
    if((class(dim(dat))!="integer")|(class(dat)[1]=="data.frame")){
      isWarning = TRUE
      messages = c(messages,"\nWARNING: counts are not provided in AIT.anndata$raw[['X']]. Counts are **REQUIRED** for the schema, but most downstream scrattch.taxonomy and scrattch.mapping functions should still work.")
    } else {
      messages = c(messages,":-) AIT.anndata$raw[['X']] looks correct.")
    }
  }
  
  ## Check the sample metadata (obs)
  required.schema.columns = c("brain_region","species")  #  UPDATE THIS LIST AS THE SCHEMA UPDATES
  if(class(AIT.anndata$obs)[1]!="data.frame"){
    isValid = FALSE
    messages = c(messages,"\nERROR: AIT.anndata$obs is not a valid data frame.")
  } else if (sum(is.element(c("cluster","cluster_label"), colnames(AIT.anndata$obs)))==0){
    isValid = FALSE
    messages = c(messages,"\nERROR: AIT.anndata$obs does not contain a column for clusters.")
  } else {
    messages = c(messages,":-) AIT.anndata$obs looks valid (additional warnings, if any, will be listed below).")
    missing.schema.columns = setdiff(required.schema.columns,gsub("_label","", colnames(AIT.anndata$obs)))
    if (length(missing.schema.columns)>0){
      isWarning = TRUE
      val = paste0(missing.schema.columns,collapse=", ")
      messages = c(messages,paste0("\nWARNING: the following AIT.anndata$obs columns are **REQUIRED** for the schema, but downstream scrattch.taxonomy and scrattch.mapping functions should still work: ",val,"."))
    }
  }
  
  ## Check the gene metadata (var)  # NOTE: THIS WILL LIKELY NEED TO BE UPDATED
  if(class(AIT.anndata$var)[1]!="data.frame"){
    isValid = FALSE
    messages = c(messages,"\nERROR: AIT.anndata$var is not a valid data frame.")
  } else if (sum(is.element(c("highly_variable_genes"), colnames(AIT.anndata$var)))==0){  # Revisit if this can be a warning instead of an Error
    isValid = FALSE
    messages = c(messages,"\nERROR: AIT.anndata$var does not contain highly_variable_genes, which is required for generating UMAPs and dendrograms.")
  } else {
    messages = c(messages,":-) AIT.anndata$var looks valid (additional warnings, if any, will be listed below).")
    if (length(grep("nsembl", colnames(AIT.anndata$var)))==0){
      isWarning = TRUE
      messages = c(messages,"\nWARNING: Ensembl IDs for genes are **REQUIRED** for the schema in AIT.anndata$var. This will not impact scrattch.mapping functionality, but will negatively impact interactions with cellxgene and other external tools.")
    }
    if (length(grep("arker", colnames(AIT.anndata$var)))==0){  # This may need to be updated later
      messages = c(messages,"Currently marker genes are not being stored in AIT.anndata$var. Storing marker genes here is not required.")
    }
  }
  
  ## Check for a 2D UMAP / latent space (obsm)  
  if((class(dim(AIT.anndata$obsm$umap))!="integer")|(class(AIT.anndata$obsm$umap)[1]=="data.frame")|(!is.element(class(AIT.anndata$obsm$umap[1,1]),c("integer","numeric","character")))){
    isWarning = TRUE
    messages = c(messages,"\nWARNING: A UMAP is invalid or not provided in AIT.anndata$obsm$umap.")
  }
  
  #########################################
  ## Check the "uns" sections for each mode
  
  ## Check general uns categories
  required.schema.columns = c("title")  #  UPDATE THIS LIST AS THE SCHEMA UPDATES
  missing.schema.columns = setdiff(required.schema.columns,names(AIT.anndata$uns))
  if (length(missing.schema.columns)>0){
    isWarning = TRUE
    val = paste0(missing.schema.columns,collapse=", ")
    messages = c(messages,paste0("\nWARNING: the following AIT.anndata$uns variables are **REQUIRED** for the schema, but downstream scrattch.taxonomy and scrattch.mapping functions should still work: ",val,"."))
  }
  
  ## NOTE: we may want to remove this warning section
  required.schema.columns = c("clusterInfo","clustersUse") 
  missing.schema.columns = setdiff(required.schema.columns,names(AIT.anndata$uns))
  if (length(missing.schema.columns)>0){
    isWarning = TRUE
    val = paste0(missing.schema.columns,collapse=", ")
    messages = c(messages,paste0("\nWARNING: the following AIT.anndata$uns variables are **REQUIRED** for the downstream scrattch.taxonomy and scrattch.mapping functions, but should have been generated with buildTaxonomy(): ",val,"."))
  } # End potential removal section
  
  ## Check taxonomy directory (AIT.anndata$uns$taxonomyDir)
  dat <- AIT.anndata$uns$taxonomyDir
  if(class(dat)!="character"){
    isValid = FALSE
    messages = c(messages,"\nERROR: A directory location must exist as a character file path in AIT.anndata$uns$taxonomyDir.")
  } else if(!file.exists(file.path(dat))) {
    isValid = FALSE
    messages = c(messages,"\nWARNING: the folder",dat,"is not found.")
  } else if (substr(dat,1,1)=="\\"){
    isWarning = TRUE
    messages = c(messages,"\nWARNING: the folder",dat,"should have a UNIX file structure not a Windows file structure.")
  } else {
    messages = c(messages,paste(":-) AIT.anndata$uns$taxonomyDir looks correct:",dat))
  }
  
  ## Now check the modes
  modes <- names(AIT.anndata$uns$filter)
  if (length(modes)==0){
    isValid = FALSE
    messages = c(messages,"\nERROR: taxonomy modes with filters are not found. Scrattch.taxonomy requires at least a standard mode with all cells included.  Likely this h5ad is an earlier version of scrattch.taxonomy format and should be remade. No additional parts of AIT.anndata$uns will be tested for validation.")
  } else {
    if(!is.element(AIT.anndata$uns$mode, modes)){
        isWarning = TRUE
        messages = c(messages,"\nWARNING: The current taxonomy mode is not one of the modes with available filters. Run mappingMode to change.")
    }
    for (mode in modes){
      messages = c(messages,paste("\n====== Reviewing AIT.anndata$uns for mode",mode,"=====."))
      
      ## Check the dendrogram is correct and can be loaded in from json format
      tryCatch({
          dat = json_to_dend(AIT.anndata$uns$dend[[AIT.anndata$uns$mode]])
        },
        error=function(cond) {
          isValid = FALSE
          messages = c(messages,"\nERROR: the dendrogram", dat, "is not correct please rebuild taxonomy or check with creator.")
        },
        finally={
          messages = c(messages,paste0(":-) AIT.anndata$uns$filter looks correct for mode ",mode,"."))
      })    
      
      ## Check the filter AIT.anndata$uns$filter[[mode]]
      dat <- AIT.anndata$uns$filter[[mode]]
      if(is.null(dat)){
        isValid = FALSE
        messages = c(messages,"\nERROR: the AIT.anndata$uns$filter for mode",mode,"does not exist.")
      } else if(class(dat[1])!="logical"){
        isWarning = TRUE
        messages = c(messages,"\nWARNING: the AIT.anndata$uns$filter for mode",mode,"is not a logical vector. This will likely cause problems for other scrattch functions.")  # Maybe needs to be an error
      } else if (sum(!dat)==0){
        isWarning = TRUE
        messages = c(messages,"\nWARNING: the filter for mode",mode,"excludes all of the data!")
      } else {
        messages = c(messages,paste0(":-) AIT.anndata$uns$filter looks correct for mode ",mode,"."))
      }
      
      ## Check the QC_markers AIT.anndata$uns$filter[[mode]]
      dat <- AIT.anndata$uns$QC_markers[[mode]]
      if(is.null(dat)){
        messages = c(messages,paste0("No QC_markers are calculated for mode ",mode,"."))
      } else{
        required.inputs = c("allMarkers", "classBr", "countsQC", "cpmQC", "markers", "qc_genes", "qc_samples", "subclassF")  
        missing.inputs = setdiff(required.inputs,names(dat))
        if (length(missing.inputs)>0){
          val = paste0(missing.inputs,collapse=", ")
          messages = c(messages,paste0("\nThe following AIT.anndata$uns$QC_markers variables are missing: ",val,". IF THIS MODE IS INTENDED FOR TREE MAPPING, please run buildPatchseqTaxonomy() function with mode set as ",mode," to calculate these values."))
        } else {
          messages = c(messages,paste0(":-) QC_markers are calculated and look correct for mode ",mode,"."))
        }
      }
    }
  }
 
  ## Output/return results
  write(messages,file.path(log.file.path,"checkTaxonomy_log.txt"))
  if(!isValid) {
    warning("AIT.anndata is not in scrattch.taxonomy format. See checkTaxonomy_log.txt for details.")
  } else if(isWarning) {
    warning("AIT.anndata is in scrattch.taxonomy format, BUT some warnings were seen. See checkTaxonomy_log.txt for details.")
  } else {
    print("AIT.anndata appears to be a completely valid scrattch.taxonomy file. You should be good to go!")
  }
  return(isValid)
}

