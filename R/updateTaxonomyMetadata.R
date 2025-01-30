#' Checks whether an anndata object is in scrattch.taxonomy format and returns a log-file if not
#'
#' @param metadata A metadata table (data.frame) to be included in the obs slot of an AIT file that will follow the AIT schema
#' @param log.file.path The directory to output the logfile of errors and warnings (if any; default getwd())
#' @param standardize_metadata If TRUE (default) will clean up standard schema files to try and remove common errors (e.g., differences in case, trailing spaces, etc.), and converting to factors.  
#' @param compute_ontology_terms A vector with any of the following terms: "ontology_term_id" for "organism", "anatomical_region", "self_reported_sex", "self_reported_ethnicity", "assay", or "disease". By default (NULL) no terms will be computed. For any terms included, will look for that column name and will return the best fit ontology ids in a new column with "ontology_term_id" appended.
#' @param compute_brain_atlas_terms Default (NULL) skips this step. If provided, can be one of: DHBA (developing human brain atlas), HBA (human brain atlas), or MBA (mouse brain atlas), which correspond to ontologies of the same name at https://github.com/brain-bican. 
#' @param compute_CL_terms Default (NULL) skips this step and is strongly recommended.  If a column name is provided (e.g., "subclass") will attempt to find the nearest CL term for whatever if included in that column. This will only provide reasonable results if this column includes human readable names that are similar to names found in cell ontology.
#' @param compute_ensembl_terms If TRUE (not recommended) will attempt to convert from gene symbol to Ensembl ID. Default is FALSE.
#' 
#' @import stringdist
#' @import ontologyIndex
#' 
#' Note: Any breaking issues will cause this function to return FALSE.  And potential issues will still return TRUE but will output a warning to stderr.  All messages will get returned to the log file. 
#' 
#' @return A list where "metadata" is the updated metadata file and there are addition list entries corresponding to any additional columns and statistics around the conversions performed in compute_ontology_terms, if any
#'
#' @export
updateTaxonomyMetadata = function(metadata, 
                                  log.file.path = getwd(),
                                  standardize_metadata = TRUE,
                                  compute_brain_atlas_terms = NULL,
                                  compute_ontology_terms = NULL,
                                  compute_CL_terms = NULL,
                                  compute_ensembl_terms = FALSE)
  {
  #########################################
  ## Initial general check and set up
  
  if(sum(grepl("data.frame",tolower(class(metadata))))==0){
    stop("metadata is not a variable of class data.frame.")
  }
  
  messages   = NULL
  output     = list(starting_metadata=metadata)
  
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
    ## For this section, we'll just go through one variable at a time for things that we'd want to adjust
    
    ## cluster_id and [cellannotation_setname]
    # NO CHECKS NOW
    # In future, could do some CAS lookup table, but probably would be a separate function
    
    ## load_id, donor_id, assay
    # NO CHECKS NOW
    # Should be standard values provided from upstream BKP processing
    
    column.name = "organism"
    # Should be standard values provided from upstream BKP processing
    # For now converting to sentence case since this is how organism scientific names are stored in NCBI
    if(column.name %in% schema.columns){
      # data(species) # this is in case we need a way to convert from common ("human") to scientific ("homo sapiens") names
      metadata[,column.name] <- tosentence(metadata[,column.name])  # will convert back to factor later
    }
    
    ## donor_age
    # Need to follow up on this. Is it numeric or categorical?
    
    column.name = "anatomic_region"
    # Convert the first character to lowercase - this is the way the vast majority of UBERON terms look and will make it easier to look up ontology terms
    if(column.name %in% schema.columns){
      metadata[,column.name] <- firsttolower(metadata[,column.name])  # will convert back to factor later
    }
    
    ## self_reported_sex
    # NO CHECKS NOW
    # Should be standard values provided from upstream BKP processing
    
    ## self_reported_ethnicity
    # NO CHECKS NOW
    # Might be standard values provided from upstream BKP processing, but also likely lot's of blanks.
    
    ## disease
    # Should add a check to convert anything likely healthy, control, normal, etc. to a single term (e.g., "control"?)
    
    column.name = "suspension_type"
    # First convert to lower case, then check if the distance of the values to "cell" or "nucleus" are less than 2 (to account for misspellings) and set everything else to "na"
    if(column.name %in% schema.columns){
      value    <-  as.character(metadata[,column.name])
      value    <- tolower(value)
      celldist <- stringdist("cell", value, method = "lv")
      nucdist  <- stringdist("nucleus", value, method = "lv")
      value[celldist<=2]     <- "cell"
      value[nucdist<=2]      <- "nucleus"
      value[!is.element(value,c("cell","nucleus"))] = "na"
      metadata[,column.name] <- value  # will convert back to factor later
    }
  
    column.name = "is_primary_data"
    # Check if boolean (logical) if not, convert to one by looking for terms similar to false and setting the rest to true
    if(column.name %in% schema.columns){
      if(is.logical(metadata[,column.name])){
        # No nothing if it's a logical... this is what we want!
      } else {
      value     <-  as.character(metadata[,column.name])
      value     <- tolower(value)
      falsedist <- stringdist("false", value, method = "lv")
      value[falsedist<=2]    <- "false"
      metadata[,column.name] <- (value!="false") # if "false" set to FALSE, otherwise set to TRUE
    }
    
  }

  #########################################
  ## If desired, compute ontology terms

  ## Start with terms that have a single obo-formatted ontology to work with  
  ontology_terms = c("organism", "assay", "anatomical_region", "self_reported_ethnicity", "disease", "cell_type")
  ontology_prefix= setNames(c("NCBITaxon:","EFO","UBERON","HANCESTRO","MONDO","CL"),ontology_terms)
  ontology_urls  = setNames(c("https://raw.githubusercontent.com/obophenotype/ncbitaxon/refs/heads/master/subsets/taxslim.obo",
                              "http://www.ebi.ac.uk/efo/efo.obo",
                              "https://github.com/obophenotype/uberon/releases/download/v2022-08-19/amniote-basic.obo",
                              "https://raw.githubusercontent.com/EBISPOT/hancestro/main/hancestro.obo",
                              "https://github.com/monarch-initiative/mondo/releases/download/v2022-09-06/mondo.obo",
                              "https://purl.obolibrary.org/obo/cl.obo"),
                            ontology_terms)
  ontology_files = setNames(c("ncbitaxon.obo","efo.obo","uberon.obo","hancestro.obo","mondo.obo","cl.obo"),ontology_terms)
  
  ontology_terms_to_compute = setdiff(interect(ontology_terms,compute_ontology_terms),"cell_type")  # CL is done separately below
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
        ontology <- ontologyIndex::get_OBO(ontology_files[onto_term])
        kp       <- substr(ontology$id,1,nchar(ontology_prefix[onto_term]))==ontology_prefix[onto_term]
        for (i in 1:length(ontology)) 
          ontology[[i]] <- ontology[[i]][kp]
        ## Find the best match ontology id for the relevant metadata column that contains names
        ontology_vector     <- as.character(metadata[,onto_term])
        unique_onto_terms   <- unique(ontology_vector)  # We only need to look up each term once
        ontology_conversion <- as.data.frame(data.table::rbindlist(lapply(unique_onto_terms,._find_best_ontology_match,ontology)))
        ontology_conversion <- ontology_conversion[match(ontology_vector,unique_onto_terms)]
        ontology_conversion$original_name <- ontology_vector
        ## Save the returned ID, returned name, string distance from original name, and original name for each item
        output[[onto_term]] = ontology_conversion
        ## Add a new column to the metadata table for the ontology ID terms
        metadata[,paste0(onto_term,"_ontology_term_id")] = ontology_conversion$id
        ## Report any changed names
        if(sum(ontology_conversion$distance)>0){
          change <- paste0(ontology_conversion[,1],"('",ontology_conversion[,4],"'==>'",ontology_conversion[,2],"')")
          messages = c(messages, paste0("\nWARNING: the following ontology terms do not correspond perfectly with inputted values:\n"))
          messages = c(messages, paste0("-----",paste(sort(unique(change[ontology_conversion$distance>0])),collapse="\n-----")))
        }
      }
    }
  }
  
  ## Now deal with "self_reported_sex"
  if ("self_reported_sex" %in% compute_ontology_terms){
    ontology_vector <- as.character(metadata[,"self_reported_sex"])
    out_vector <- rep("unknown", length(ontology_vector))
    out_vector[tolower(ontology_vector)=="female"]  = "PATO_0000383"
    out_vector[tolower(ontology_vector)=="male"]    = "PATO_0000384"
    out_vector[tolower(ontology_vector)=="f"]       = "PATO_0000383" # Deal with common abbreviations
    out_vector[tolower(ontology_vector)=="m"]       = "PATO_0000384" # Deal with common abbreviations
    metadata[,"self_reported_sex_ontology_term_id"] = out_vector
  }
  
  
  ## Now deal with "brain_region"
  
  ### NOTE: I might need a table converting between abbrevations and longer names if we are using abbreviations.  Right now this only looks at the longer names.  Return synonym
  
  # Available ontologies listed here
  brain_atlas_terms = c("DHBA", "HBA", "MBA")
  brain_atlas_prefix= setNames(c("DHBA", "HBA", "MBA"),brain_atlas_terms)
  brain_atlas_urls  = setNames(c("https://github.com/brain-bican/developing_human_brain_atlas_ontology/raw/refs/heads/main/dhbao-simple-non-classified.obo",
                                 "https://github.com/brain-bican/human_brain_atlas_ontology/raw/refs/heads/main/hbao-simple-non-classified.obo",
                                 "https://github.com/brain-bican/mouse_brain_atlas_ontology/raw/refs/heads/main/mbao-simple-non-classified.obo"),
                               brain_atlas_terms)
  brain_atlas_files = setNames(c("dhbao.obo","hbao.obo","mbao.obo"),brain_atlas_terms)
  skip_message = FALSE
  if (is.null(compute_brain_atlas_terms)) skip_message = TRUE
  compute_brain_atlas_terms = intersect(compute_brain_atlas_terms,brain_atlas_terms)

  ## Now do the comparison
  if (is.null(compute_brain_atlas_terms)){
    if(!skip_message)
      messages = c(messages, paste0("\nWARNING: compute_brain_atlas_terms are not included in existing options. Skipping brain mapping."))
  } else {
    ## Search for existing downloaded obo file and use that if in current directory, otherwise download
    onto_term = compute_brain_atlas_terms[1]
    if(!file.exists(brain_atlas_files[onto_term])){
      file <- try(download.file(brain_atlas_urls[onto_term],brain_atlas_files[onto_term])) 
      if("try-error" %in% class(file))
        messages = c(messages, paste0("\nWARNING: ",onto_term," ontology file not available and is being skipped"))
    } 
    if(file.exists(brain_atlas_files[onto_term])){
      ## Read in the ontology file
      ontology <- ontologyIndex::get_OBO(brain_atlas_files[onto_term]) # Need to add extract_tags="everything") to get synonym, 
      # Need to reformat synonyms to search them.  Example: [1] "\"NP\" EXACT []"
      kp       <- substr(ontology$id,1,nchar(brain_atlas_prefix[onto_term]))==brain_atlas_prefix[onto_term]
      for (i in 1:length(ontology)) 
        ontology[[i]] <- ontology[[i]][kp]
      ## Find the best match ontology id for the relevant metadata column that contains names
      ontology_vector     <- as.character(metadata[,onto_term])
      unique_onto_terms   <- unique(ontology_vector)  # We only need to look up each term once
      ontology_conversion <- as.data.frame(data.table::rbindlist(lapply(unique_onto_terms,._find_best_ontology_match,ontology)))
      ontology_conversion <- ontology_conversion[match(ontology_vector,unique_onto_terms),]
      ontology_conversion$original_name <- ontology_vector
      ## Save the returned ID, returned name, string distance from original name, and original name for each item
      output[["brain_region"]] = ontology_conversion
      ## Add a new column to the metadata table for the ontology ID terms
      metadata[,"brain_region_ontology_term_id"] = ontology_conversion$id
      ## Report any changed names
      if(sum(ontology_conversion$distance)>0){
        change <- paste0(ontology_conversion[,1],"('",ontology_conversion[,4],"'==>'",ontology_conversion[,2],"')")
        messages = c(messages, paste0("\nWARNING: the following ontology terms do not correspond perfectly with inputted values:\n"))
        messages = c(messages, paste0("-----",paste(sort(unique(change[ontology_conversion$distance>0])),collapse="\n-----")))
      }
    }
  }
  
  
  ## Now deal with "cell_type_ontology_term"
  
  if(!is.null(compute_CL_terms)){
    if (compute_CL_terms[1] %in% schema.columns){
      onto_term <- "cell_type"
      ## Search for existing downloaded obo file and use that if in current directory, otherwise download
      if(!file.exists(ontology_files[onto_term])){
        file <- try(download.file(ontology_urls[onto_term],ontology_files[onto_term])) 
        if("try-error" %in% class(file))
          messages = c(messages, paste0("\nWARNING: ",onto_term," ontology file not available and is being skipped"))
      } 
      if(file.exists(ontology_files[onto_term])){
        ## Read in the ontology file
        ontology <- ontologyIndex::get_OBO(ontology_files[onto_term])
        kp       <- substr(ontology$id,1,nchar(ontology_prefix[onto_term]))==ontology_prefix[onto_term]
        for (i in 1:length(ontology)) 
          ontology[[i]] <- ontology[[i]][kp]
        ## Find the best match cell ontology id for the relevant metadata column in compute_CL_terms
        ontology_vector     <- as.character(metadata[,compute_CL_terms[1]]) 
        unique_onto_terms   <- unique(ontology_vector)  # We only need to look up each term once
        ontology_conversion <- as.data.frame(data.table::rbindlist(lapply(unique_onto_terms,._find_best_ontology_match,ontology)))
        ontology_conversion <- ontology_conversion[match(ontology_vector,unique_onto_terms)]
        ontology_conversion$original_name <- ontology_vector
        ## Save the returned ID, returned name, string distance from original name, and original name for each item
        output[[onto_term]] = ontology_conversion
        ## Add a new column to the metadata table for the ontology ID terms
        metadata[,"cell_type_ontology_term"] = ontology_conversion$id
        ## Report any changed names
        if(sum(ontology_conversion$distance)>0){
          change <- paste0(ontology_conversion[,1],"('",ontology_conversion[,4],"'==>'",ontology_conversion[,2],"')")
          messages = c(messages, paste0("\nWARNING: the following ontology terms do not correspond perfectly with inputted values:\n"))
          messages = c(messages, paste0("-----",paste(sort(unique(change[ontology_conversion$distance>0])),collapse="\n-----")))
        }
      }
    } else {
      messages = c(messages, paste0("\nWARNING: column name ",compute_CL_terms[1]," is not found so cell_type_ontology_term compute is being skipped."))
    }
  }
  
  
  #########################################
  ## Convert any categorical variables to factors, if they aren't already
  
  ## NOTE: This part needs to be updated to retain the original factor information, if available. I'll create a new helper function for this that will take two variables (orignal factor vector and new character vector) as input.
  
  
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







#' Convert to sentence case
#'
#' @param str A character string
#'
#' @return The character string in sentence case
tosentence <- function(str) {
  str <- as.character(str)
  paste0(toupper(substr(str,1,1)),tolower(substr(str,2,nchar(str))))
}


#' Convert the first character of a string to lowercase
#'
#' @param str A character string
#'
#' @return The character string with specific characters changed to lower case
firsttolower <- function(str) {
  str <- as.character(str)
  paste0(tolower(substr(str,1,1)),substr(str,2,nchar(str)))
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