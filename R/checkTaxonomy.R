#' Checks whether an anndata object is in scrattch.taxonomy format and returns a log-file if not
#'
#' Any breaking issues will cause this function to return FALSE.  And potential issues will still return TRUE but will output a warning to stderr.  All messages will get returned to the log file. 
#'
#' @param AIT.anndata A reference taxonomy anndata object to be tested
#' @param log.file.path The directory to output the logfile of errors and warnings (if any; default getwd())
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

  data(schema)
  
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
    missing.schema.columns = c()
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
    if(is.element(element, names(AIT.anndata$uns))){
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
    validation = ._validate_schema_element(AIT.anndata$uns[[element]], column_def, messages, isValid)
    messages = validation[["messages"]]; isValid = validation[["isValid"]]
  }

  ## Now we will check any RECOMMENDED schema elements that are present in uns
  recommended.schema.columns = ._get_schema_elements(schema, "uns", "RECOMMENDED")
  tovalidate.schema.columns = intersect(recommended.schema.columns, colnames(AIT.anndata$uns)) ## FIX or something very close to this.
  for(element in tovalidate.schema.columns){
    if(is.element(elemenet, names(AIT.anndata$uns))){
      column_def = ._get_schema_def(element)
      validation = ._validate_schema_element(AIT.anndata$uns[[element]], column_def, messages, isValid)
      messages = validation[["messages"]]; isValid = validation[["isValid"]]
    }
  }

  #########################################
  ## Check the metadata (var)
  ## highly_variable_genes[_name], marker_genes[_name], ensembl_id all of which are RECOMMENDED

  validation = .validate_var_elements(AIT.anndata, messages, isValid, isWarning)
  messages = validation[["messages"]]; isValid = validation[["isValid"]]; isWarning = validation[["isWarning"]]

  #########################################
  ## Check the embeddings X_[embedding] which are RECOMMENDED

  validation = ._validated_embeddings(AIT.anndata, messages, isValid, isWarning)
  messages = validation[["messages"]]; isValid = validation[["isValid"]]; isWarning = validation[["isWarning"]]

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
    ## @Jeremy Should we print the messages to the console as well?
  } else if(isWarning) {
    warning("AIT.anndata is in Allen Institute Taxonomy format, BUT some warnings were seen. See checkTaxonomy_log.txt for details.")
  } else {
    print("AIT.anndata appears to be a completely valid Allen Institute Taxonomy file!")
  }

  ## Return the logical for validation
  return(isValid)
}

#' This function will return information about a given schema Key
#'
#' @param key The schema Key to provide information about to the user.
#'
#' @return internal schema information
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
#' @return schema elements
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
#' @param column The vector of entrees in the column to validate.
#' @param column_def The schema definition for the column.
#' @param messages The current messages to append to.
#' @param isValid The current call for whether the taxonomy is valid.
#' @param pull_cl If FALSE (default) loads a preset list of CL terms from OBO; otherwise pulls from OBO
#' @param validate_percent_cl Percent of entries that must correspond to valid CL terms to validate
#' @param pull_assay If FALSE (default) loads the list of EFO terms (assays); otherwise pulls from EBI (VERY slow) 
#' @param pull_ncbitaxon If FALSE (default) loads the list of species with gene information at NCBI; otherwise pulls from OBO (VERY slow) 
#' @param pull_uberon If FALSE (default) loads the list of anatomic regions from UBERON; otherwise pulls from OBO
#' @param pull_brain_atlases If FALSE (default) loads the list of brain atlas ids; otherwise pulls from brain-bican
#' @param pull_hancestro If FALSE (default) loads the list of HANCESTRO terms; otherwise, pulls from OBO 
#' @param pull_mondo If FALSE (default) loads the list of MONDO terms; otherwise, pulls from OBO 
#' @param pull_ensembl If FALSE (default) loads the list of Ensembl terms from NCBI; otherwise, pulls from NCBI (VERY slow)
#' @param validate_percent_ensembl Percent of entries that must correspond to valid Ensembl terms to validate
#'
#' @return A list with messages and an isValid logical call
#'
#' @keywords internal
._validate_schema_element = function(column, column_def, messages, isValid,
                                     pull_assay = FALSE,
                                     pull_cl = FALSE, validate_percent_cl = 80,            # cell_type_ontology_term variables
                                     pull_ncbitaxon = FALSE,                               # organism_ontology_term_id variables
                                     pull_uberon = FALSE,                                  # anatomical_region_ontology_term_id variables
                                     pull_brain_atlases = FALSE,                           # brain_region_ontology_term_id variables
                                     pull_hancestro = FALSE,                               # self_reported_ethnicity_ontology_term_id variables
                                     pull_mondo = FALSE,                                   # disease_ontology_term_id variables
                                     pull_ensembl = FALSE, validate_percent_ensembl = 60   # ensembl_id variables
                                     ){

  ###########################################################
  ## First we will validate individual columns type matches expectation from schema

  if(column_def$Type == "str"){
      if(!all(is.character(column))){
        messages = c(messages, paste0("\nERROR: The anndata.obs element: ", column_def$Key, " must be strings."))
        isValid = FALSE
      }
  }else if(column_def$Type == "bool"){
      if(!all(column %in% c(TRUE, FALSE))){ 
        messages = c(messages, paste0("\nERROR: The anndata.obs element: ", column_def$Key, " must be boolean."))
        isValid = FALSE
      }
  }else if(column_def$Type == "Categorical"){
      if(!all(is.factor(column))){
        messages = c(messages, paste0("\nERROR: The anndata.obs element: ", column_def$Key, " must be categorical (factor)."))
        isValid = FALSE
      }
  }

  ###########################################################
  ## Now we will validate individual columns against schema

  ## suspension_type
  if(column_def$Key == "suspension_type"){
      if(!all(column %in% c("cell", "nucleus", "na"))){
        messages = c(messages, paste0("\nERROR: The anndata.obs element: ", column_def$Key, " must be either 'cell', 'nucleus' or 'na'."))
        isValid = FALSE
      }
  }

  
  ## cell_type_ontology_term
  ## Must exist as a CL term to be valid. If unknown use 'CL:0000003' for native cell.
  data(cl)
  if(pull_cl){
    # This step requires an internet connection and could take several minutes
    file   <- try(download.file("https://purl.obolibrary.org/obo/cl.obo","cl.obo"))
    if("try-error" %in% class(file)){
      messages = c(messages, paste0("\nWARNING: CL terms not accessible via https://purl.obolibrary.org/obo/cl.obo and not updated."))
    } else {
      cl_obo <- ontologyIndex::get_OBO("cl.obo") 
      cl     <- cl_obo$id[substr(cl_obo$id,1,2)=="CL"]
    }
  }
  # Allow for NA values
  cl <- c(cl,"NA")
  column[is.na(column)] = "NA"
  # Now do the test
  if(column_def$Key == "cell_type_ontology_term"){
    if(!all(column %in% cl)){
      percent_with_cl  <- signif(100*mean(column %in% cl),4)
      message  = paste0("\nThe anndata.obs element: ", column_def$Key, " contains ",percent_with_cl,"% CL terms.")
      if(percent_with_cl<validate_percent_cl){
        message = paste0(message,"\nERROR: At least ",validate_percent_cl,"% CL terms required to validate.")
        isValid = FALSE
      }
      messages = c(messages, paste0(message))
    }
  }
  
  
  ## organism_ontology_term_id
  #  Must be a child of http://purl.obolibrary.org/obo/NCBITaxon_33208 (for Metazoa); e.g., "NCBITaxon:9606" for Homo sapiens
  #  Default file includes all 658 species with gene info at NCBI
  data(ncbitaxon)
  if(pull_ncbitaxon){
    #file   <- try(download.file("http://purl.obolibrary.org/obo/ncbitaxon.obo","ncbitaxon.obo")) # complete, HUGE file
    file   <- try(download.file("https://raw.githubusercontent.com/obophenotype/ncbitaxon/refs/heads/master/subsets/taxslim.obo","ncbitaxon.obo"))
    if("try-error" %in% class(file)){
      messages = c(messages, paste0("\nWARNING: NCBITaxon terms not accessible via https://purl.obolibrary.org/obo/ncbitaxon.obo and not updated."))
    } else {
      ncbitaxon_obo <- ontologyIndex::get_OBO("ncbitaxon.obo")
      ncbitaxon     <- unique(c(ncbitaxon,ncbitaxon_obo$id[substr(ncbitaxon_obo$id,1,10)=="NCBITaxon:"])) # 
    }
  }
  # Allow for NA values
  ncbitaxon <- c(ncbitaxon,"NA")
  column[is.na(column)] = "NA"
  # Now do the test
  if(column_def$Key == "organism_ontology_term_id"){
    if(!all(column %in% ncbitaxon)){
      messages = c(messages, paste0("\nERROR: The anndata.obs element: ", column_def$Key, " must all be valid NCBITaxons."))
      isValid = FALSE
    }
  }
 
   
  ## assay_ontology_term_id
  #  Most appropriate EFO ontology term for assay. (e.g., "10x 3' v2"="EFO:0009899", "10x 3' v3"="EFO:0009922", "Smart-seq"="EFO:0008930")
  #  Note: using the latest release as of 1/15/2025 
  data(assay)
  if(pull_assay){
    file   <- try(download.file("http://www.ebi.ac.uk/efo/efo.obo","efo.obo"))
    if("try-error" %in% class(file)){
      messages = c(messages, paste0("\nWARNING: NCBITaxon terms not accessible via https://github.com/EBISPOT/efo/ and not updated."))
    } else {
      efo_obo <- ontologyIndex::get_OBO("efo.obo")
      assay   <- efo_obo$id[substr(efo_obo$id,1,3)=="EFO"]
    }
  }
  # Allow for NA values
  assay <- c(assay,"NA")
  column[is.na(column)] = "NA"
  # Now do the test
  if(column_def$Key == "assay_ontology_term_id"){
    if(!all(column %in% assay)){
      messages = c(messages, paste0("\nERROR: The anndata.obs element: ", column_def$Key, " must all be valid EFO terms."))
      isValid = FALSE
    }
  }
  
  
  ## anatomical_region_ontology_term_id
  # Must be an UBERON term (http://www.obofoundry.org/ontology/uberon.html) to be valid <- NOTE: may need to update to allow for Allen Brain Map ontology term
  # NOTE: This code will need to be updated for more recent UBERON terms, but the current one matches CELLxGENE
  data(uberon)
  if(pull_uberon){
    file   <- try(download.file("https://github.com/obophenotype/uberon/releases/download/v2022-08-19/amniote-basic.obo","uberon.obo"))
      # NOTE: we are using and older version of UBERON because the one at http://purl.obolibrary.org/obo/uberon.obo has a cycle error.
    if("try-error" %in% class(file)){
      messages = c(messages, paste0("\nWARNING: UBERON terms not accessible via https://purl.obolibrary.org/obo/uberon.obo and not updated."))
    } else {
      uberon_obo <- ontologyIndex::get_OBO("uberon.obo")
      uberon     <- uberon_obo$id[substr(uberon_obo$id,1,6)=="UBERON"]
    }
  }
  # Allow for NA values
  uberon <- c(uberon,"NA")
  column[is.na(column)] = "NA"
  # Now do the test
  if(column_def$Key == "anatomical_region_ontology_term_id"){
    if(!all(column %in% uberon)){
      messages = c(messages, paste0("\nERROR: The anndata.obs element: ", column_def$Key, " must all be valid NCBITaxons."))
      isValid = FALSE
    }
  }
  
  
  ## brain_region_ontology_term_id
  # If provided, must be a valid id from developing_human_brain_atlas_ontology, human_brain_atlas_ontology, or mouse_brain_atlas_ontology
  # These are the Allen Institute / BICAN ontologies available on brain-bican, and could be expanded in the future
  data(brain_atlases)
  if(pull_brain_atlases){
    # Developing human brain atlas
    file <- try(download.file("https://github.com/brain-bican/developing_human_brain_atlas_ontology/raw/refs/heads/main/dhbao-simple-non-classified.obo","dhbao.obo"))
    if("try-error" %in% class(file)){
      messages = c(messages, paste0("\nWARNING: DHBA terms not accessible via https://github.com/brain-bican/ and not updated."))
    } else {
      dhbao_obo <- ontologyIndex::get_OBO("dhbao.obo")
      brain_atlases <- c(brain_atlases, dhbao_obo$id[substr(dhbao_obo$id,1,4)=="DHBA"])
    }
    # Human brain atlas
    file <- try(download.file("https://github.com/brain-bican/human_brain_atlas_ontology/raw/refs/heads/main/hbao-simple-non-classified.obo","hbao.obo"))
    if("try-error" %in% class(file)){
      messages = c(messages, paste0("\nWARNING: HBA terms not accessible via https://github.com/brain-bican/ and not updated."))
    } else {
      hbao_obo <- ontologyIndex::get_OBO("hbao.obo")
      brain_atlases <- c(brain_atlases, hbao_obo$id[substr(hbao_obo$id,1,3)=="HBA"])
    }
    # Mouse brain atlas
    file <- try(download.file("https://github.com/brain-bican/mouse_brain_atlas_ontology/raw/refs/heads/main/mbao-simple-non-classified.obo","mbao.obo"))
    if("try-error" %in% class(file)){
      messages = c(messages, paste0("\nWARNING: MBA terms not accessible via https://github.com/brain-bican/ and not updated."))
    } else {
      mbao_obo <- ontologyIndex::get_OBO("mbao.obo")
      brain_atlases <- c(brain_atlases, mbao_obo$id[substr(mbao_obo$id,1,3)=="MBA"])
    }
  }
  # Allow for NA values
  brain_atlases <- c(brain_atlases,"NA")
  column[is.na(column)] = "NA"
  # Now do the test
  if(column_def$Key == "brain_region_ontology_term_id"){
    if(!all(column %in% brain_atlases)){
      messages = c(messages, paste0("\nERROR: The anndata.obs element: ", column_def$Key, " must all be valid brain atlas terms from BICAN-brain ontologies."))
      isValid = FALSE
    }
  }
  
  
  ## self_reported_sex_ontology_term_id
  # female = PATO_0000383; male = PATO_0000384; other uncommon options are children of http://purl.obolibrary.org/obo/PATO_0001894
  # for now, we set to NA anything that doesn't fit into one of the above two categories (but it will be retained in self_reported_sex) <-- THIS IS NOT IDEAL
  # Called "sex_ontology_term_id" in cellxgene
  # For now we validate anything, but report the number of rows that don't have PATO_0000383 or PATO_0000384. Again, not ideal, but ok for now.
  if(column_def$Key == "self_reported_sex_ontology_term_id"){
    column[is.na(column)] = "NA"
    if(!all(column %in% c("PATO_0000383", "PATO_0000384"))){
      percent_MF  <- signif(100*mean(column %in% c("PATO_0000383", "PATO_0000384")),4)
      messages = c(messages, paste0("\nWARNING: The anndata.obs element: ", column_def$Key, " contains ",percent_MF,"% of terms as PATO_0000383 or PATO_0000384. This may cause issues in translation to CELLxGENE."))
    }
  }
  
  
  ## self_reported_ethnicity_ontology_term_id
  # If organism_ontolology_term_id is "NCBITaxon:9606" for Homo sapiens, this MUST be either a HANCESTRO term, "multiethnic" if more than one ethnicity is reported, or "unknown" if unavailable. 
  # For now we have a parameter asking if the species is "NCBITaxon:9606" and only flag as invalid if TRUE and above requirements not met
  # We are looking in a slightly different place for this ontology
  data(hancestro)
  if(pull_hancestro){
    file   <- try(download.file("https://raw.githubusercontent.com/EBISPOT/hancestro/main/hancestro.obo","hancestro.obo"))
    if("try-error" %in% class(file)){
      messages = c(messages, paste0("\nWARNING: HANCESTRO terms not accessible via https://raw.githubusercontent.com/EBISPOT/hancestro/main/hancestro.obo and not updated."))
    } else {
      hancestro_obo <- ontologyIndex::get_OBO("hancestro.obo")
      hancestro     <- hancestro_obo$id[substr(hancestro_obo$id,1,9)=="HANCESTRO"]
    }
  }
  hancestro <- c(hancestro,"multiethnic","unknown")
  # Allow for NA values
  hancestro <- c(hancestro,"NA")
  column[is.na(column)] = "NA"
  # Now do the test
  if(column_def$Key == "organism_ontology_term_id"){
    if(!all(column %in% hancestro)){
      messages = c(messages, paste0("\nERROR: The anndata.obs element: ", column_def$Key, " must all be valid HANCESTRO terms, or 'multiethnic' or 'unknown'."))
      isValid = FALSE
    }
  }
  
  
  ## disease_ontology_term_id
  # Either "PATO_0000461" or a "MONDO" term
  # Note that http://purl.obolibrary.org/obo/mondo.obo has a cycle so we are backtracking to the version used in cellxgene.  Should update eventually.
  data(mondo)
  if(pull_mondo){
    file   <- try(download.file("https://github.com/monarch-initiative/mondo/releases/download/v2022-09-06/mondo.obo","mondo.obo")) 
    if("try-error" %in% class(file)){
      messages = c(messages, paste0("\nWARNING: HANCESTRO terms not accessible via https://github.com/monarch-initiative/mondo/releases/download/v2022-09-06/mondo.obo and not updated."))
    } else {
      mondo_obo <- ontologyIndex::get_OBO("mondo.obo")
      mondo     <- mondo_obo$id[substr(mondo_obo$id,1,5)=="MONDO"]
    }
  }
  mondo <- c(mondo,"PATO:0000461")
  # Allow for NA values
  mondo <- c(mondo,"NA")
  column[is.na(column)] = "NA"
  # Now do the test
  if(column_def$Key == "organism_ontology_term_id"){
    if(!all(column %in% mondo)){
      messages = c(messages, paste0("\nERROR: The anndata.obs element: ", column_def$Key, " must all be valid HANCESTRO terms, or 'multiethnic' or 'unknown'."))
      isValid = FALSE
    }
  }
  
  
  ## cluster_algorithm
  # How should we check this? We need to write a schema on https://github.com/AllenInstitute/AllenInstituteTaxonomy/tree/main/schema

  ## 
  return(list(messages=messages, isValid=isValid))
}

#' This function will validate the modes in an AIT file
#'
#' @param AIT.anndata The AIT anndata object to validate.
#' @param messages The current messages to append to.
#' @param isValid The current isValid status.
#' @param isWarning The current isWarning status.
#'
#' @return A list with messages, an isValid logical call, and an isWarning logical call
#'
#' @keywords internal
._validate_ait_modes = function(AIT.anndata, messages, isValid, isWarning){
  
  ## Gather all modes
  modes <- names(AIT.anndata$uns$filter)

  ## Check modes
  if (length(modes)==0){
    isValid = FALSE
    messages = c(messages,"\nERROR: taxonomy modes with filters are not found. Allen Institute Taxonomy requires at least a standard mode with all cells included which should have been created with `buildTaxonomy`.  Likely this h5ad is an earlier version of Allen Institute Taxonomy (AIT) format and should be remade.")
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
          messages = c(messages,paste0("\n:-) AIT.anndata$uns$dend looks correct for mode ", mode, "."))
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
        messages = c(messages, paste0("\n:-) AIT.anndata$uns$filter looks correct for mode ",mode,"."))
      }
    }
  }
  return(list("messages"= messages, "isValid" = isValid, "isWarning" = isWarning))
}

#' This function will validate the embeddings in an AIT file
#'
#' @param AIT.anndata The AIT anndata object to validate.
#' @param messages The current messages to append to.
#' @param isValid The current isValid status.
#' @param isWarning The current isWarning status.
#'
#' @return A list with messages, an isValid logical call, and an isWarning logical call
#'
#' @keywords internal
._validated_embeddings = function(AIT.anndata, 
                                   messages, 
                                   isValid, 
                                   isWarning){

  ## Check for a 2D UMAP / latent space (obsm)
  embeddings = names(AIT.anndata$obsm)

  if(length(embeddings) > 0){
    for(embedding in embeddings){
      
      ## Check the naming scheme for embeddings which should follow "X_[embedding]"
      if(!grepl("^X_", embedding)){
        isValid = FALSE
        messages = c(messages,"\nERROR: An embedding is provided but not in the correct format of 'X_[embedding]'")
      } else {
        messages = c(messages,"\n:-) AIT.anndata$obsm contains correctly named embeddings (additional warnings, if any, will be listed below).")
      }

      ## Now check that each embedding is what we would expect
      if((class(dim(AIT.anndata$obsm[[embedding]])) != "integer") | 
          (class(AIT.anndata$obsm[[embedding]])[1]=="data.frame") |
          (!is.element(class(AIT.anndata$obsm[[embedding]][1,1]), c("integer","numeric","character")))){
        isWarning = TRUE
        messages = c(messages, paste0("\nWARNING: An embedding: ", embedding, " is invalid."))
      }
    }
  }

  return(list("messages"= messages, "isValid" = isValid, "isWarning" = isWarning))
}

#' This function will validate the var in an AIT file
#'
#' @param AIT.anndata The AIT anndata object to validate.
#' @param messages The current messages to append to.
#' @param isValid The current isValid status.
#' @param isWarning The current isWarning status.
#'
#' @return A list with messages, an isValid logical call, and an isWarning logical call
#'
#' @keywords internal
.validate_var_elements = function(AIT.anndata, messages, isValid, isWarning, validate_percent_ensembl=60, pull_ensembl=FALSE){

  ## Check the highly_variable_genes
  if(sum(is.element(c("highly_variable_genes"), colnames(AIT.anndata$var)))==0){
    isWarning = TRUE
    messages = c(messages,"\nWARNING: AIT.anndata$var does not contain highly_variable_genes[_name], which is recommended for generating UMAPs and dendrograms.")
  }else{
    var_gene_columns = colnames(AIT.anndata$var)[grepl("highly_variable_genes", colnames(AIT.anndata$var))]
    for(var_gene_column in var_gene_columns){
      column = AIT.anndata$var[[var_gene_column]]
      if(!all(column %in% c(TRUE, FALSE))){ 
        isValid = FALSE
        messages = c(messages, paste0("\nThe anndata.var element: ", var_gene_column, " must be boolean."))
      }else{
        messages = c(messages,"\n:-) AIT.anndata$var contains highly_variable_genes (additional warnings, if any, will be listed below).")
      }
    }
  }

  ## Check the marker_genes
  if(sum(is.element(c("marker_genes"), colnames(AIT.anndata$var)))==0){
    isWarning = TRUE
    messages = c(messages,"\nWARNING: AIT.anndata$var does not contain marker_genes[_name], which is recommended for generating UMAPs and dendrograms.")
  }else{
    marker_gene_columns = colnames(AIT.anndata$var)[grepl("marker_genes", colnames(AIT.anndata$var))]
    for(marker_gene_column in marker_gene_columns){
      column = AIT.anndata$var[[marker_gene_column]]
      if(!all(column %in% c(TRUE, FALSE))){ 
        isValid = FALSE
        messages = c(messages, paste0("The anndata.var element: ", marker_gene_column, " must be boolean.\n"))
      }else{
        messages = c(messages,"\n:-) AIT.anndata$var contains marker_genes (additional warnings, if any, will be listed below).")
      }
    }
  }
  
  
  ## Check the ensembl_id
  if(is.element("ensembl_id", colnames(AIT.anndata$var))){
    isWarning = TRUE
    messages = c(messages,"\nWARNING: AIT.anndata$var does not contain ensembl_id, which is recommended for generating UMAPs and dendrograms.")
  }else{
    ## Validate ensembl_id???
    messages = c(messages,"\n:-) AIT.anndata$var contains ensembl_id (additional warnings, if any, will be listed below).")
  
    ## At least 60% of terms must exist to be valid. 
    ## Will look at predefined list for "human", "mouse", "marmoset", and "rhesus_macaque" by default and will only query other species if asked.
    data(ensembl)
    if(pull_ensembl){
      ncbi_gene_info <- try(data.table::fread("https://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2ensembl.gz"))
      if("try-error" %in% class(ncbi_gene_info)){
        messages = c(messages, paste0("\nWARNING: ensembl terms not accessible via https://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2ensembl.gz and not updated.\n"))
      } else {
        ensembl  <- convertEns$Ensembl_gene_identifier
      }
    }

    ## Allow for NA values
    ensembl <- c(ensembl, "NA")

    ## Pull ensembl_id column out
    column.to.validate = AIT.anndata$var$ensembl_id
    column.to.validate[is.na(column.to.validate)] = "NA"

    ## Now do the test
    if(!all(column.to.validate %in% ensembl)){
      percent_with_ensembl  <- signif(100*mean(column.to.validate %in% ensembl),4)
      message  = paste0("\nThe anndata.obs element: ", "ensembl_id", " contains ", percent_with_ensembl,"% ensembl terms.")
      if(percent_with_ensembl < validate_percent_ensembl){
        message = paste0(message,"\n:ERROR: At least ",validate_percent_ensembl,"% ensembl terms required to validate.")
        isValid = FALSE
      }
      messages = c(messages, paste0(message))
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