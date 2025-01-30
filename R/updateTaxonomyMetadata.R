## Compute most likely "ontology_term_id" for "organism", "anatomical_region", "brain_region", "self_reported_sex", "self_reported_ethnicity", "assay", "disease", and "ensembl_id"
## Also deal with "na"/"NA"/NA and other similar issues

#' Checks whether an anndata object is in scrattch.taxonomy format and returns a log-file if not
#'
#' @param metadata A metadata table (data.frame) to be included in the obs slot of an AIT file that will follow the AIT schema
#' @param log.file.path The directory to output the logfile of errors and warnings (if any; default getwd())
#' @param standardize_metadata If TRUE (default) will clean up standard schema files to try and remove common errors (e.g., differences in case, trailing spaces, etc.), and converting to factors.  
#' @param compute_ontology_terms A vector with any of the following terms: "ontology_term_id" for "organism", "anatomical_region", "brain_region", "self_reported_sex", "self_reported_ethnicity", "assay", or "disease". By default (NULL) no terms will be computed. For any terms included, will look for that column name and will return the best fit ontology ids in a new column with "ontology_term_id" appended.
#' @param compute_ensembl If TRUE (not recommended) will attempt to convert from gene symbol to Ensembl ID. Default is FALSE.
#' 
#' Note: Any breaking issues will cause this function to return FALSE.  And potential issues will still return TRUE but will output a warning to stderr.  All messages will get returned to the log file. 
#' 
#' @return A list where "metadata" is the updated metadata file and there are addition list entries corresponding to any additional columns and statistics around the conversions performed in compute_ontology_terms, if any
#'
#' @export
updateTaxonomyMetadata = function(metadata, 
                                  log.file.path = getwd(),
                                  standardize_metadata = TRUE,
                                  compute_ontology_terms = NULL,
                                  compute_ensembl = FALSE)
  {
  #########################################
  ## Initial general check and set up
  
  if(sum(grepl("data.frame",tolower(class(metadata))))==0){
    stop("metadata is not a variable of class data.frame.")
  }
  
  messages   = NULL
  output     = list()
  
  ## Read in the schema columns for consideration in this function
  data(schema)
  all.schema.columns = c(._get_schema_elements(schema, "obs"),._get_schema_elements(schema, "obs", "RECOMMENDED"))
  schema.columns = intersect(all.schema.columns, colnames(metadata))
  
  ## Check that desired ontology terms to compute exist in the table
  missing_ontology_terms = setdiff(compute_ontology_terms, schema.columns) 
  if(length(missing_ontology_terms)>0)
    messages = c(messages, paste0("\nWARNING: ",paste(missing_ontology_terms,collapse=", ")," are missing from either the schema or the metadata table and cannot be computed."))
  compute_ontology_terms = intersect(compute_ontology_terms, schema.columns) 
  
  
  #########################################
  ## If desired, standardize the metadata
  if(standardize_metadata){
    # For this section, we'll just go through one variable at a time for things that we'd want to adjust
    
    # Start with organism
    if("organism" %in% schema.columns){
      ## THIS IS STILL NOT COMPLETED, AND MAY NOT BE NECESSARY
      data(species) # this is a way to convert from common ("human") to scientific ("homo sapiens") names
    }
  
  }

  #########################################
  ## If desired, compute ontology terms

  ## Start with terms that have a single obo-formatted ontology to work with  
  ontology_terms = c("organism", "anatomical_region", "self_reported_ethnicity", "assay", "disease")
  ontology_urls  = setNames(c("https://raw.githubusercontent.com/obophenotype/ncbitaxon/refs/heads/master/subsets/taxslim.obo",
                              "http://www.ebi.ac.uk/efo/efo.obo",
                              "https://github.com/obophenotype/uberon/releases/download/v2022-08-19/amniote-basic.obo",
                              "https://raw.githubusercontent.com/EBISPOT/hancestro/main/hancestro.obo",
                              "https://github.com/monarch-initiative/mondo/releases/download/v2022-09-06/mondo.obo"),
                            ontology_terms)
  ontology_files = setNames(c("ncbitaxon.obo","efo.obo","uberon.obo","hancestro.obo","mondo"),ontology_terms)
  
  ontology_terms_to_compute = interect(ontology_terms,compute_ontology_terms)
  if(length(ontology_terms_to_compute)>0){
    for (onto_term in ontology_terms_to_compute){
      ## Search for existing downloaded obo file and use that if in current directory, otherwise download
      if(!file.exists(ontology_files[onto_term])){
        file <- try(download.file(ontology_urls[onto_term],ontology_files[onto_term])) 
        if("try-error" %in% class(file))
          messages = c(messages, paste0("\nWARNING: ",onto_term," ontology file not available and is being skipped"))
      } 
      if(file.exists(ontology_files[onto_term])){
        ## Read in the ontology file
        ontology            <- ontologyIndex::get_OBO(ontology_files[onto_term])
        ## Find the best match ontology id for the relevant metadata column that contains names
        ontology_vector     <- as.character(metadata[,onto_term])
        ontology_conversion <- as.data.frame(data.table::rbindlist(lapply(ontology_vector,._find_best_ontology_match,ontology)))
        ontology_conversion$original_name <- ontology_vector
        ## Save the returned ID, returned name, string distance from original name, and original name for each item
        output[[onto_term]] = ontology_conversion
        ## Add a new column to the metadata table for the ontology ID terms
        metadata[,paste0(onto_term,"_ontology_term_id")] = ontology_conversion$id
      }
    }
  }
  
  ## Now deal with "self_reported_sex"
  if ("self_reported_sex" %in% compute_ontology_terms){
    ontology_vector <- as.character(metadata[,"self_reported_sex"])
    out_vector <- rep("unknown", length(ontology_vector))
    out_vector[tolower(ontology_vector)=="female"]  = "PATO_0000383"
    out_vector[tolower(ontology_vector)=="male"]    = "PATO_0000384"
    out_vector[tolower(ontology_vector)=="f"]       = "PATO_0000383" # Deal with abbreviations
    out_vector[tolower(ontology_vector)=="m"]       = "PATO_0000384" # Deal with abbreviations
    metadata[,"self_reported_sex_ontology_term_id"] = out_vector
  }
  
  
  ## Now deal with "brain_region"
  # LEAVING THIS BLANK RIGHT NOW.  
  if ("brain_region" %in% compute_ontology_terms){
    messages = c(messages, paste0("\nWARNING: automated brain_region prediction is not yet implemented and is being skipped"))
  }
  # Currently the data in the files below are not formatted in a way that I can access the brain region names
  # "https://github.com/brain-bican/mouse_brain_atlas_ontology/raw/refs/heads/main/mbao-simple-non-classified.obo"
  # "https://github.com/brain-bican/human_brain_atlas_ontology/raw/refs/heads/main/hbao-simple-non-classified.obo"
  # "https://github.com/brain-bican/developing_human_brain_atlas_ontology/raw/refs/heads/main/dhbao-simple-non-classified.obo"
  # There are larger, more complete files but they are REALLY slow (e.g., I have not been successful in even opening them with R)
  
  
  #########################################
  ## Convert any categorical variables to factors, if they aren't already
  
  schema.columns = intersect(all.schema.columns, colnames(metadata))
  for(element in schema.columns){
    column_def = ._get_schema_def(element)
    if(column_def$Type == "Categorical")
      if(!("factor" %in% class(metadata[,element])))
        metadata[,element] <- as.factor(metadata[,element])
  }
  
  
  #########################################
  ## Write out the log file and return the updated metadata table
  
  ## Output/return results
  if(!is.null(messages)){
    write(messages,file.path(log.file.path,"checkMetadata_log.txt"))
    print("\n=== Messages from updateTaxonomyMetadata.R are available in checkMetadata_log.txt. ===\n")
    ## @Jeremy Should we print the messages to the console as well?
  } 

  ## Return the metadata information
  output[["metadata"]] = metadata
  return(output)
  
}





#' This function will return the best matching ontology term id given a name and a distance to that best match
#'
#' @param search_term Name of an ontology term you want to try and match (e.g., "frontal cortex" or "human")
#' @param ontology Variable of class ontology_index that holds a complete ontology read in from obo format using get_OBO.
#' @param top_n Number of top matches to return (default = 1)
#'
#' @import stringdist
#' @import ontologyIndex
#'
#' @return character indicating the associated id name and distance for the closest matching term
#'
#' @keywords internal
._find_best_ontology_match = function(search_term, ontology, top_n=1){
  
  # Extract term names and IDs
  term_names <- ontology$name
  term_ids <- names(term_names)  # IDs are stored as names of the name vector
  
  # Calculate string similarity (e.g., Levenshtein distance)
  distances <- stringdist(search_term, term_names, method = "lv")
  
  # Find the top N matches (lowest distances)
  top_indices <- order(distances)[1:min(top_n, length(term_names))]
  
  # Return the matching IDs and names
  matches <- data.frame(
    id = term_ids[top_indices],
    name = term_names[top_indices],
    distance = distances[top_indices]
  )
  
  return(matches)
  
}