#' Updates a metadata data frame with ontology terms
#'
#' See updateTaxonomyMetadata for details. computeOntologyTerms is a wrapper for updateTaxonomyMetadata with defaults to do compute the ontology terms for everything (except CL), but not do anything else.
#' Any breaking issues will cause this function to return FALSE.  And potential issues will still return TRUE but will output a warning to stderr.  All messages will get returned to the log file. 
#'
#' @import stringdist
#' @import ontologyIndex
#' 
#' @return A list where "metadata" is the updated metadata file and there are additional list entries corresponding to conversions surrounding ontology terms.  The original metadata data frame is also returned.
#'
#' @export
computeOntologyTerms <- function(
    metadata, 
    log.file.path             = getwd(),
    log.file.name             = "computeOntologyTerms_log.txt",
    standardize.metadata      = FALSE,
    compute.ontology.terms    = c("organism", "anatomical_region", "self_reported_sex", "self_reported_ethnicity", "assay", "disease"),
    compute.brain.atlas.terms = "DHBA",
    convert.regions.to.names  = TRUE,
    compute.cl.terms          = NULL,
    print.messages            = FALSE){
  ## Call updateTaxonomyMetadata immediately
  updateTaxonomyMetadata(metadata                  = metadata,
                         log.file.path             = log.file.path,
                         log.file.name             = log.file.name,
                         compute.ontology.terms    = compute.ontology.terms,
                         compute.brain.atlas.terms = compute.brain.atlas.terms,
                         convert.regions.to.names  = convert.regions.to.names,
                         compute.cl.terms          = compute.cl.terms,
                         print.messages            = print.messages)
}


#' Updates a metadata data frame to better align with the AIT schema
#'
#' This function defaults to standardizing metadata but not messing with ontologies. computeOntologyTerms is a wrapper function that defaults to compute the ontology terms for everything (except CL), but not do anything else.  Either function can be used identically by adjusting the parameters.
#' Any breaking issues will cause this function to return FALSE.  And potential issues will still return TRUE but will output a warning to stderr.  All messages will get returned to the log file. 
#'
#' @param metadata A metadata table (data.frame) to be included in the obs slot of an AIT file that will follow the AIT schema
#' @param log.file.path The directory to output the logfile of errors and warnings (if any; default getwd())
#' @param standardize.metadata If TRUE (default) will clean up standard schema files to try and remove common errors (e.g., differences in case, trailing spaces, etc.), and converting to factors.  
#' @param compute.ontology.terms A vector with any of the following terms: "ontology_term_id" for "organism", "anatomical_region", "self_reported_sex", "self_reported_ethnicity", "assay", or "disease". By default (NULL) no terms will be computed. For any terms included, will look for that column name and will return the best fit ontology ids in a new column with "ontology_term_id" appended.
#' @param compute.brain.atlas.terms Default (NULL) skips this step. If provided, can be one of: DHBA (developing human brain atlas), HBA (human brain atlas), or MBA (mouse brain atlas), which correspond to ontologies of the same name at https://github.com/brain-bican. 
#' @param convert.regions.to.names If full brain region names are provided in "anatomical_region", this does nothing. Otherwise, updateTaxonomyMetadata will attempt to convert brain region abbreviations to brain region names, as required to convert to UBERON and brain atlas ontologies. If TRUE, these brain region names will overwrite inputted brain region abbrevations in "anatomical_region"; otherwise these are only returned in list entries for brain region-related ontology terms (default). Note that this variable requires a value for compute.brain.atlas.terms (not NULL), as that is the ontology that it will use to try and convert between abbreviation and name.
#' @param compute.cl.terms Default (NULL) skips this step and is strongly recommended.  If a column name is provided (e.g., "subclass") will attempt to find the nearest CL term for whatever if included in that column. This will only provide reasonable results if this column includes human readable names that are similar to names found in cell ontology.
#' @param print.messages Print messages only to a log file (FALSE; default) or also to the screen (TRUE)
#' 
#' @import stringdist
#' @import ontologyIndex
#' 
#' @return A list where "metadata" is the updated metadata file and there are addition list entries corresponding to any additional columns and statistics around the conversions performed in compute.ontology.terms, if any
#'
#' @export
updateTaxonomyMetadata = function(metadata, 
                                  log.file.path             = getwd(),
                                  log.file.name             = "updateTaxonomyMetadata_log.txt",
                                  standardize.metadata      = TRUE,
                                  compute.ontology.terms    = NULL,
                                  compute.brain.atlas.terms = NULL,
                                  convert.regions.to.names  = FALSE,
                                  compute.cl.terms          = NULL,
                                  print.messages            = FALSE)
  {
  
  #########################################
  ## Initial general check and set up
  
  if(sum(grepl("data.frame",tolower(class(metadata))))==0){
    stop("metadata is not a variable of class data.frame.")
  }
  
  ## Set up messages and output variables (where the initial metadata file is stored)
  messages   = NULL
  output     = list(starting_metadata=metadata)
  
  ## List of ontologies that have a single obo-formatted ontology to work with  
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
  
  # Now list the available brain region atlas ontologies listed here
  brain_atlas_terms = c("DHBA", "HBA", "MBA")
  brain_atlas_prefix= setNames(c("DHBA", "HBA", "MBA"),brain_atlas_terms)
  brain_atlas_urls  = setNames(c("https://github.com/brain-bican/developing_human_brain_atlas_ontology/raw/refs/heads/main/dhbao-simple-non-classified.obo",
                                 "https://github.com/brain-bican/human_brain_atlas_ontology/raw/refs/heads/main/hbao-simple-non-classified.obo",
                                 "https://github.com/brain-bican/mouse_brain_atlas_ontology/raw/refs/heads/main/mbao-simple-non-classified.obo"),
                               brain_atlas_terms)
  brain_atlas_files = setNames(c("dhbao.obo","hbao.obo","mbao.obo"),brain_atlas_terms)
  
  ## Read in the schema columns for consideration in this function
  data(schema)
  all.schema.columns = c(._get_schema_elements(schema, "obs"),._get_schema_elements(schema, "obs", "RECOMMENDED"))
  schema.columns = intersect(all.schema.columns, colnames(metadata))
  
  ## Check that desired ontology terms to compute exist in the table
  missing_ontology_terms = setdiff(compute.ontology.terms, schema.columns) 
  if(length(missing_ontology_terms)>0)
    messages = c(messages, paste0("\nWARNING: ",paste(missing_ontology_terms,collapse=", ")," are missing from either the schema or the metadata table and cannot be computed."))
  compute.ontology.terms = intersect(compute.ontology.terms, schema.columns) 
 
  skip_message = FALSE  # This refers to a message below in the "brain_region_ontology_term_id" section
  if (is.null(compute.brain.atlas.terms)) skip_message = TRUE
  compute.brain.atlas.terms = intersect(compute.brain.atlas.terms,brain_atlas_terms)
  
  
  #########################################
  ## If "anatomical_region" exists, check whether terms are abbreviations and if so, convert to full brain region 'name' 
  # Assume abbreviations are provided if average of <10 characters per entry
  
  if("anatomical_region" %in% compute.ontology.terms){
    num.characters = mean(as.numeric(lapply(as.character(metadata[,"anatomical_region"]),nchar)),na.rm=TRUE)
    # frac.lowercase = mean(as.numeric(lapply(as.character(metadata[,"anatomical_region"]),function(x) stringdist(x,tolower(x))/nchar(x))),na.rm=TRUE) # Not used since several abbreviations in DHBA are all lowercase.
    is.region.abbreviations = FALSE
    if(num.characters<10) is.region.abbreviations = TRUE
    
    # If regions are abbreviations, convert to full names (we will convert BACK later if convert.regions.to.names=FALSE)
    if(is.region.abbreviations){
      if (is.null(compute.brain.atlas.terms)){
        messages = c(messages, paste0("\nWARNING: compute.brain.atlas.terms are not included in existing options, so brain region names cannot be computed from brain region abbreviations and therefore ontology definitions based on brain regions will likely be WRONG."))
      } else {
        ## Message that we are going to look for long names
        messages = c(messages, paste0("\nMESSAGE: Brain region abbreviations appear to be included in anatomical_region. Attempting to convert to full name for future ontology conversions, if requested."))
        ## Search for existing downloaded obo file and use that if in current directory, otherwise download
        onto_term = compute.brain.atlas.terms[1]
        if(!file.exists(brain_atlas_files[onto_term])){
          file <- try(download.file(brain_atlas_urls[onto_term],brain_atlas_files[onto_term])) 
          if("try-error" %in% class(file))
            messages = c(messages, paste0("\nWARNING: ",onto_term," ontology file not available and is being skipped"))
        } 
        if(file.exists(brain_atlas_files[onto_term])){
          ## Read in the ontology file and format synonyms (e.g., "\"NP\" EXACT []" --> "NP")
          ontology <- ontologyIndex::get_OBO(brain_atlas_files[onto_term],extract_tags="everything") 
          kp       <- substr(ontology$id,1,nchar(brain_atlas_prefix[onto_term]))==brain_atlas_prefix[onto_term]
          for (i in 1:length(ontology)) 
            ontology[[i]] <- ontology[[i]][kp]
          synonym  <- as.character(lapply(ontology$synonym,function(x) x[1]))  # To deal with 0 or 2+ synonyms
          synonym  <- substr(gsub("\"","",synonym),1,nchar(synonym)-11)        # Formatting DHBA/HBA/MBA
          ontology$name <- setNames(synonym, as.character(ontology$name))      # Formatting for use with ._find_best_ontology_match
          
          ## Find the best match ontology name for the relevant metadata column that contains abbreviations (synonyms)
          ontology_vector     <- as.character(metadata[,"anatomical_region"])
          unique_onto_terms   <- unique(ontology_vector)  # We only need to look up each term once
          ontology_conversion <- as.data.frame(data.table::rbindlist(lapply(unique_onto_terms,._find_best_ontology_match,ontology)))
          ontology_conversion <- ontology_conversion[match(ontology_vector,unique_onto_terms),]
          ontology_conversion$original_name <- ontology_vector
          colnames(ontology_conversion) <- c("name","synonym","distance","original_synonym")
          ## Save the returned ID, returned name, string distance from original name, and original name for each item
          output[["brain_region_abbreviation_conversion"]] = ontology_conversion
          ## Add a new column to the metadata table for the ontology ID terms
          metadata[,"anatomical_region"] = ontology_conversion$name
          ## Report any changed names
          if(sum(ontology_conversion$distance)>0){
            change <- paste0(ontology_conversion[,1],"('",ontology_conversion[,4],"'==>'",ontology_conversion[,2],"')")
            messages = c(messages, paste0("\nWARNING: the following ontology synonyms do not correspond perfectly with inputted values:\n"))
            messages = c(messages, paste0("-----",paste(sort(unique(change[ontology_conversion$distance>0])),collapse="\n-----")))
          }
        }
      }
    }  
  } #else {
    #is.region.abbreviations = FALSE
    #anatomical_region_name = metadata[,"anatomical_region"]
  #}
  
  
  #########################################
  ## If desired, standardize the metadata
  
  if(standardize.metadata){
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
    
    column.name = "anatomical_region"
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
      # Do nothing if it's a logical... this is what we want!
      if(!is.logical(metadata[,column.name])){
        value     <-  as.character(metadata[,column.name])
        value     <- tolower(value)
        falsedist <- stringdist("false", value, method = "lv")
        value[falsedist<=2]    <- "false"
        metadata[,column.name] <- (value!="false") # if "false" set to FALSE, otherwise set to TRUE
      }
    }
  }

  #########################################
  ## If desired, compute ontology terms
    
  ## Start with all the default terms
  ontology_terms_to_compute = setdiff(intersect(ontology_terms,compute.ontology.terms),"cell_type")  # CL is done separately below
  if(length(ontology_terms_to_compute)>0){
    for (onto_term in ontology_terms_to_compute){
      if(onto_term %in% colnames(metadata)){
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
          ontology_conversion <- ontology_conversion[match(ontology_vector,unique_onto_terms),]
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
      } else {
        messages = c(messages, paste0("\nWARNING: ",onto_term," not include in metadata so ontology terms not calculated."))
      }
    }
  }
  
    
  ## Now deal with "self_reported_sex"
  if ("self_reported_sex" %in% compute.ontology.terms){
    ontology_vector <- as.character(metadata[,"self_reported_sex"])
    out_vector <- rep("unknown", length(ontology_vector))
    out_vector[tolower(ontology_vector)=="female"]  = "PATO_0000383"
    out_vector[tolower(ontology_vector)=="male"]    = "PATO_0000384"
    out_vector[tolower(ontology_vector)=="f"]       = "PATO_0000383" # Deal with common abbreviations
    out_vector[tolower(ontology_vector)=="m"]       = "PATO_0000384" # Deal with common abbreviations
    metadata[,"self_reported_sex_ontology_term_id"] = out_vector
  }
  
  
  ## Now deal with "brain_region"
  if (is.null(compute.brain.atlas.terms)|(!("anatomical_region" %in% schema.columns))){
    if(!skip_message)
      messages = c(messages, paste0("\nWARNING: compute.brain.atlas.terms are not included in existing options. Skipping brain mapping."))
  } else {
    ## Search for existing downloaded obo file and use that if in current directory, otherwise download
    onto_term = compute.brain.atlas.terms[1]
    if(!file.exists(brain_atlas_files[onto_term])){
      file <- try(download.file(brain_atlas_urls[onto_term],brain_atlas_files[onto_term])) 
      if("try-error" %in% class(file))
        messages = c(messages, paste0("\nWARNING: ",onto_term," ontology file not available and is being skipped"))
    } 
    if(file.exists(brain_atlas_files[onto_term])){
      ## Read in the ontology file
      ontology <- ontologyIndex::get_OBO(brain_atlas_files[onto_term]) 
      # Need to reformat synonyms to search them.  Example: [1] "\"NP\" EXACT []"
      kp       <- substr(ontology$id,1,nchar(brain_atlas_prefix[onto_term]))==brain_atlas_prefix[onto_term]
      for (i in 1:length(ontology)) 
        ontology[[i]] <- ontology[[i]][kp]
      ## Find the best match ontology id for the relevant metadata column that contains names
      ontology_vector     <- as.character(metadata[,"anatomical_region"])
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
  if(!is.null(compute.cl.terms)){
    if (compute.cl.terms[1] %in% schema.columns){
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
        ## Find the best match cell ontology id for the relevant metadata column in compute.cl.terms
        ontology_vector     <- as.character(metadata[,compute.cl.terms[1]]) 
        unique_onto_terms   <- unique(ontology_vector)  # We only need to look up each term once
        ontology_conversion <- as.data.frame(data.table::rbindlist(lapply(unique_onto_terms,._find_best_ontology_match,ontology)))
        ontology_conversion <- ontology_conversion[match(ontology_vector,unique_onto_terms),]
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
      messages = c(messages, paste0("\nWARNING: column name ",compute.cl.terms[1]," is not found so cell_type_ontology_term compute is being skipped."))
    }
  }
  
  
  #########################################
  ## Convert any categorical variables to factors, if they aren't already
  
  ## If convert.regions.to.names=FALSE, and "anatomical_region" exists change "anatomical_region" back to original values
  initial.metadata <- output[["starting_metadata"]]
  if((!convert.regions.to.names)&("anatomical_region" %in% schema.columns)&("anatomical_region" %in% colnames(metadata))){
    metadata[,"anatomical_region"] <- initial.metadata[,"anatomical_region"]
  }
    
  ## Factorize initial columns that are part of the schema
  initial.schema.columns  <- intersect(all.schema.columns, colnames(initial.metadata))
  if(length(initial.schema.columns)>0) 
    for(element in initial.schema.columns){
      column_def = ._get_schema_def(element)
      if(column_def$Type == "Categorical")
        metadata[,element] <- ._transfer_factor_levels(metadata[,element], initial.metadata[,element])
  }
  
  ## Factorize other initial columns not part of the schema that were originally factors
  initial.nonschema.columns <- setdiff(colnames(initial.metadata),initial.schema.columns)
  if(length(initial.nonschema.columns)>0) 
    for(element in initial.nonschema.columns)
      if(is.factor(initial.metadata[,element]))
        metadata[,element] <- ._transfer_factor_levels(metadata[,element], initial.metadata[,element])
  
  ## Factorize ontology columns added by this function
  initial.nonschema.columns <- setdiff(colnames(initial.metadata),initial.schema.columns)
  if(length(compute.ontology.terms)>0) 
    for(element in compute.ontology.terms){
      element.id <- paste0(element,"_ontology_term_id")
      metadata[,element.id] <- ._transfer_factor_levels(metadata[,element.id],metadata[,element])
    }
  if("cell_type_ontology_term" %in% colnames(metadata))
    if(!is.null(compute.cl.terms)){
      element    <- compute.cl.terms[1]
      element.id <- "cell_type_ontology_term"
      metadata[,element.id] <- ._transfer_factor_levels(metadata[,element.id],metadata[,element])
    }
  if("brain_region_ontology_term_id" %in% colnames(metadata)){
    element    <- "anatomical_region"
    element.id <- "brain_region_ontology_term_id"
    metadata[,element.id] <- ._transfer_factor_levels(metadata[,element.id],metadata[,element])
  }


  #########################################
  ## Write out the log file and return the updated metadata table
  
  ## Output/return results
  if(!is.null(messages)){
    write(messages,file.path(log.file.path,log.file.name))
    print("=== Messages from updateTaxonomyMetadata.R are available in checkMetadata_log.txt. ===")

    ## If requested, output messages to screen
    if(print.messages) writeLines(readLines(file.path(log.file.path,log.file.name)))
    
  } else {
    print("No messages printed from updateTaxonomyMetadata.R.")
  }

  ## Return the metadata information
  output[["metadata"]] = metadata
  return(output)
  
}
######## END OF MAIN FUNCTION #########################################################

  

#' Convert one string back to a factor with matched, but different, values
#'
#' This function is speficially targeting the changes.  When a vector gets minor updates due to misspellings or ontologies, this factor order from the original vector is retained, even though the values are slightly changed.
#'
#' @param new.string A character string to convert to a factor
#' @param initial.factor A matched original factor variable to transfer to the string
#'
#' @return The character string with specific characters changed to lower case
._transfer_factor_levels <- function(new.string, initial.factor) {
  if(is.character(initial.factor))
    initial.factor <- as.factor(initial.factor)  # To work in situations where the inputs are characters
  if(is.numeric(initial.factor))
    initial.factor <- as.factor(initial.factor)  # To work in situations where the inputs are numeric (probably unlikely)
 
  level.order = match(as.character(levels(initial.factor)),as.character(initial.factor))
  new.factor  = factor(new.string,levels=unique(new.string[level.order]))
  return(new.factor)  
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
