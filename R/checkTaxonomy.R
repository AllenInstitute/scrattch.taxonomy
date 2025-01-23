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
  ## Check the expression data 
  
  ## Check the raw data in anndata.raw.X
  if(!is.null(AIT.anndata$raw[["X"]])){
    if((class(dim(AIT.anndata$raw[["X"]])) != "integer")){
      isValid = FALSE ## We have to be strict with the counts matrix, in my opinion.
      isWarning = TRUE
      messages = c(messages,"\nWARNING: It appears normalized counts are provided in AIT.anndata$raw[['X']]. Counts are **REQUIRED** for the schema, but most downstream scrattch.taxonomy and scrattch.mapping functions should still work.")
    } else {
      messages = c(messages,":-) AIT.anndata$raw[['X']] looks correct.")
    }
  }

  ## Check the normalized data in anndata.X
  if(!is.null(AIT.anndata$X)){
    if(max(AIT.anndata$X)>20) {
      isWarning = TRUE
      messages = c(messages,"\nWARNING: AIT.anndata$X has high values.  Please confirm this is log-normalized.")
    } else {
      messages = c(messages,":-) AIT.anndata$X looks correct.")
    }
  }

  #########################################
  ## Check the metadata (obs)
  
  ## Check that the sample metadata (obs) all exist.
  required.schema.columns = ._get_schema_elements(schema, "obs")
  if(sum(is.element(required.schema.columns, colnames(AIT.anndata$obs))) < length(required.schema.columns)){
    missing.schema.columns = setdiff(required.schema.columns, colnames(AIT.anndata$obs))
    if(length(missing.schema.columns) > 0){
      isValid = FALSE
      val = paste0(missing.schema.columns, collapse=", ")
      messages = c(messages,paste0("\nWARNING: the following AIT.anndata$obs columns are **REQUIRED** for the schema: ", val, "."))
    }
  }else{
    messages = c(messages,":-) AIT.anndata$obs contains all required schema elements (additional warnings, if any, will be listed below).")
  }

  ## If the user has provide some/all schema columns lets validate against the schema.
  tovalidate.schema.columns = setdiff(required.schema.columns, missing.schema.columns)
  for(element in tovalidate.schema.columns){
    column_def = ._get_schema_def(element)
    validation = ._validate_schema_element(AIT.anndata$obs[[element]], column_def, messages, isValid)
    messages = validation[["messages"]]; isValid = validation[["isValid"]]
  }

  ## Now we will check any RECOMMENDED schema elements that are present in obs
  recommended.schema.columns = ._get_schema_elements(schema, "obs", "RECOMMENDED")
  for(element in tovalidate.schema.columns){
    if(is.element(elemenet, names(AIT.anndata$uns))){
      column_def = ._get_schema_def(element)
      validation = ._validate_schema_element(AIT.anndata$obs[[element]], column_def, messages, isValid)
      messages = validation[["messages"]]; isValid = validation[["isValid"]]
    }
  }

  #########################################
  ## Check the metadata (uns)

  ## Check the unstructured metadata (uns)
  required.schema.columns = ._get_schema_elements(schema, "uns")
  if(sum(is.element(required.schema.columns, names(AIT.anndata$uns))) < length(required.schema.columns)){
    missing.schema.columns = setdiff(required.schema.columns, names(AIT.anndata$uns))
    if(length(missing.schema.columns) > 0){
      isValid = FALSE
      val = paste0(missing.schema.columns, collapse=", ")
      messages = c(messages,paste0("\nWARNING: the following AIT.anndata$uns columns are **REQUIRED** for the schema: ", val, "."))
    }
  }else{
    messages = c(messages,":-) AIT.anndata$uns contains all required schema elements (additional warnings, if any, will be listed below).")
  }

  ## If the user has provide some/all uns lets validate against the schema.
  tovalidate.schema.columns = setdiff(required.schema.columns, missing.schema.columns)
  for(element in tovalidate.schema.columns){
    column_def = ._get_schema_def(element)
    validation = ._validate_schema_element(AIT.anndata$obs[[element]], column_def, messages, isValid)
    messages = validation[["messages"]]; isValid = validation[["isValid"]]
  }

  ## Now we will check any RECOMMENDED schema elements that are present in uns
  recommended.schema.columns = ._get_schema_elements(schema, "uns", "RECOMMENDED")
  for(element in tovalidate.schema.columns){
    if(is.element(elemenet, names(AIT.anndata$uns))){
      column_def = ._get_schema_def(element)
      validation = ._validate_schema_element(AIT.anndata$obs[[element]], column_def, messages, isValid)
      messages = validation[["messages"]]; isValid = validation[["isValid"]]
    }
  }

  #########################################
  ## Check the metadata (var)

  ## Check the gene metadata (var)  # NOTE: THIS WILL LIKELY NEED TO BE UPDATED
  if(sum(is.element(c("highly_variable_genes"), colnames(AIT.anndata$var)))==0){  # Revisit if this can be a warning instead of an Error
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
  
  #########################################
  ## Check the embeddings X_[embedding]

  ## Check for a 2D UMAP / latent space (obsm)
  embeddings = names(AIT.anndata$obsm)

  if(length(embeddings) > 0){

    ## Check the naming scheme for embeddings which should follow "X_[embedding]"
    if(sum(grepl("X_",embeddings))==0){
      isValid = FALSE
      messages = c(messages,"\nERROR: An embedding is provided but not in the correct format of 'X_[embedding]'")
    } else {
      messages = c(messages,":-) AIT.anndata$obsm contains correctly named embeddings (additional warnings, if any, will be listed below).")
    }

    ## Now check that each embeddings
    if((class(dim(AIT.anndata$obsm$X_umap)) != "integer") | 
        (class(AIT.anndata$obsm$X_umap)[1]=="data.frame") |
        (!is.element(class(AIT.anndata$obsm$X_umap[1,1]), c("integer","numeric","character")))){
      isWarning = TRUE
      messages = c(messages,"\nWARNING: A UMAP is invalid or not provided in AIT.anndata$obsm$X_umap.")
    }

  }else{
    isValid = FALSE
    messages = c(messages,"\nERROR: No embeddings are provided in AIT.anndata$obsm.")
  }

  #########################################
  ## Check the modes

  validation = ._validate_ait_modes(AIT.anndata, messages, isValid, isWarning)
  messages = validation[["messages"]]; isValid = validation[["isValid"]]; isWarning = validation[["isWarning"]]
  
  #########################################
  ## Write out the log file and return isValid.
 
  ## Output/return results
  write(messages,file.path(log.file.path,"checkTaxonomy_log.txt"))
  if(!isValid) {
    warning("AIT.anndata is not in Allen Institute Taxonomy format as described at https://github.com/AllenInstitute/AllenInstituteTaxonomy/tree/main/schema. See checkTaxonomy_log.txt for details.")
    ## Should we print the messages to the console as well?
  } else if(isWarning) {
    warning("AIT.anndata is in Allen Institute Taxonomy format, BUT some warnings were seen. See checkTaxonomy_log.txt for details.")
  } else {
    print("AIT.anndata appears to be a completely valid Allen Institute Taxonomy file!")
  }

  ## Return the logical for validation
  return(isValid)
}

#' This function will return infromation about a given schema Key
#'
#' @param key The schema Key to provide information about to the user.
#'
#' @return
#'
#' @keywords internal
._get_schema_def = function(key=NULL){
    return(schema %>% filter(Key == key))
}

#' This function will return schema information for a given anndata Component and schema element Type
#'
#' @param schema The schema data.frame.
#' @param component The schema Component to provide information about.
#' @param type The schema element Type to provide information about.
#'
#' @return
#'
#' @keywords internal
._get_schema_elements = function(schema, component, type="MUST"){
    required_schema = schema %>% 
                        filter(Component == component) %>%
                        filter(Required == type) %>%
                        filter(!grepl("\\[", Key)) %>%
                        filter(Key != "index") %>%
                        pull(Key)
    return(required_schema)
}

#' This function will validate a given column against the schema
#'
#' @param column_name The name of the column to validate.
#' @param column_def The schema definition for the column.
#' @param messages The current messages to append to.
#'
#' @return
#'
#' @keywords internal
._validate_schema_element = function(column, column_def, messages, isValid){

  ###########################################################
  ## First we will validate invidual columns type matches expectation from schema

  if(column_def$Type == "str"){
      if(!all(is.character(column))){
        messages = c(messages, paste0("The anndata.obs element: ", column_def$Key, " must be strings.\n"))
        isValid = FALSE
      }
  }else if(column_def$Type == "bool"){
      if(!all(column %in% c(TRUE, FALSE))){ 
        messages = c(messages, paste0("The anndata.obs element: ", column_def$Key, " must be boolean.\n"))
        isValid = FALSE
      }
  }else if(column_def$Type == "Categorical"){
      if(all(is.factor(column))){
        messages = c(messages, paste0("The anndata.obs element: ", column_def$Key, " must be categorical (factor).\n"))
        isValid = FALSE
      }
  }

  ###########################################################
  ## Now we will validate invidual columns against schema

  ## suspension_type
  if(column_def$Key == "suspension_type"){
      if(!all(column %in% c("cell", "nucleus", "na"))){
        messages = c(messages, paste0("The anndata.obs element: ", column_def$Key, " must be either 'cell', 'nucleus' or 'na'.\n"))
        isValid = FALSE
      }
  }

  ## cell_type_ontology_term

  ## organism_ontology_term_id

  ## anatomical_region_ontology_term_id

  ## self_reported_sex_ontology_term_id

  ## disease_ontology_term_id

  ## ensembl_id

  ## cluster_[algorithm]
  ## How should we check this? We need to write a schema on https://github.com/AllenInstitute/AllenInstituteTaxonomy/tree/main/schema

  ## 
  return(list("messages" = messages, "isValid" = isValid))
}

#' This function will validate the modes in an AIT file
#'
#' @param AIT.anndata The AIT anndata object to validate.
#' @param messages The current messages to append to.
#' @param isValid The current isValid status.
#' @param isWarning The current isWarning status.
#'
#' @return
#'
#' @keywords internal
._validate_ait_modes = function(AIT.anndata, messages, isValid, isWarning){
  
  ## Gather all modes
  modes <- names(AIT.anndata$uns$filter)

  ## Check modes
  if (length(modes)==0){
    isValid = FALSE
    messages = c(messages,"\nERROR: taxonomy modes with filters are not found. Allen Institute Taxonomy requires atleast a standard mode with all cells included which should have been created with `buildTaxonomy`.  Likely this h5ad is an earlier version of Allen Institute Taxonomy (AIT) format and should be remade.")
  } else {
    if(!is.element(AIT.anndata$uns$mode, modes)){
        isWarning = TRUE
        messages = c(messages,"\nWARNING: The current taxonomy mode is not one of the modes with available filters. Run mappingMode to change.")
    }
    for (mode in modes){
      messages = c(messages,paste("\n====== Reviewing AIT.anndata$uns for mode", mode, "=====."))
      
      ## Check the dendrogram is correct and can be loaded in from json format
      tryCatch({
          dat = json_to_dend(AIT.anndata$uns$dend[[AIT.anndata$uns$mode]])
        },
        error=function(cond) {
          isValid = FALSE
          messages = c(messages,"\nERROR: the dendrogram", dat, "is not correct please rebuild taxonomy or check with creator.")
        },
        finally={
          messages = c(messages,paste0(":-) AIT.anndata$uns$dend looks correct for mode ", mode, "."))
      })    
      
      ## Check the filter AIT.anndata$uns$filter[[mode]]
      dat <- AIT.anndata$uns$filter[[mode]]
      if(is.null(dat)){
        isValid = FALSE
        messages = c(messages,"\nERROR: the AIT.anndata$uns$filter for mode",mode,"does not exist.")
      } else if(class(dat[1])!="logical"){
        isValid = FALSE
        messages = c(messages,"\nERROR: the AIT.anndata$uns$filter for mode",mode,"is not a logical vector. This will likely cause problems for other scrattch functions.")  # Maybe needs to be an error
      } else if (sum(!dat)==0){
        isWarning = TRUE
        messages = c(messages,"\nWARNING: the filter for mode", mode, "excludes all of the data!")
      } else {
        messages = c(messages, paste0(":-) AIT.anndata$uns$filter looks correct for mode ",mode,"."))
      }
    }
  }
  return(list("messages"= messages, "isValid" = isValid, "isWarning" = isWarning))
}


## Not in the schema and well handled elsewhere (loadTaxonomy)
## Check taxonomy directory (AIT.anndata$uns$taxonomyDir)
# dat <- AIT.anndata$uns$taxonomyDir
# if(class(dat)!="character"){
#   isValid = FALSE
#   messages = c(messages,"\nERROR: A directory location must exist as a character file path in AIT.anndata$uns$taxonomyDir.")
# } else if(!file.exists(file.path(dat))) {
#   isValid = FALSE
#   messages = c(messages,"\nWARNING: the folder",dat,"is not found.")
# } else if (substr(dat,1,1)=="\\"){
#   isWarning = TRUE
#   messages = c(messages,"\nWARNING: the folder",dat,"should have a UNIX file structure not a Windows file structure.")
# } else {
#   messages = c(messages,paste(":-) AIT.anndata$uns$taxonomyDir looks correct:",dat))
# }

## QC_markers has been removed from the AIT object in favor of computing on the fly what is required for patchseq.
# ## Check the QC_markers AIT.anndata$uns$filter[[mode]]
# dat <- AIT.anndata$uns$QC_markers[[mode]]
# if(is.null(dat)){
#   messages = c(messages,paste0("No QC_markers are calculated for mode ",mode,"."))
# } else{
#   required.inputs = c("allMarkers", "classBr", "countsQC", "cpmQC", "markers", "qc_genes", "qc_samples", "subclassF")  
#   missing.inputs = setdiff(required.inputs,names(dat))
#   if (length(missing.inputs)>0){
#     val = paste0(missing.inputs,collapse=", ")
#     messages = c(messages,paste0("\nThe following AIT.anndata$uns$QC_markers variables are missing: ",val,". IF THIS MODE IS INTENDED FOR TREE MAPPING, please run buildPatchseqTaxonomy() function with mode set as ",mode," to calculate these values."))
#   } else {
#     messages = c(messages,paste0(":-) QC_markers are calculated and look correct for mode ",mode,"."))
#   }
# }