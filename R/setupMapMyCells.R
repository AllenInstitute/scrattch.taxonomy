#' This function builds files needed for hierarchical mapping and stores them in the uns$hierarchical of AIT (Shiny) taxonomy.
#'
#' @param AIT_anndata A reference taxonomy anndata object.
#' @param hierarchy List of term_set_labels in the reference taxonomy ordered from most gross to most fine.
#' @param AIT_anndata_path Local file path of the AIT reference taxonomy (h5ad file).
#' @param force Boolean value indicating whether to overwrite the AIT reference taxonomy's hierarchical file for a given mode.
#' @param n_processors Number of independent worker processes to spin up.
#' @param normalization Normalization of the h5ad files; must be either 'raw' or 'log2CPM'.
#' @param tmp_dir Temporary directory for writing out the hierarchical files.
#' @param user_precomp_stats_path Alternative path to the user provided precompute stats HDF5 file. Will be generated, if not provided.
#' @param user_query_markers_path Alternative path to the user provided query markers JSON file. Will be generated, if not provided.
#' 
#' Note: this hierarchical mapping is a wrapper around cell_type_mapper,   
#'       and call's it's functions to generate needed files needed for mapping.
#' 
#' @import anndata
#'
#' @return AIT.anndata, a reference taxonomy with hierarchical files, such as precompute stats and query markers saved in uns.
#'
#' @export
addMapMyCells = function(AIT_anndata,
                             hierarchy=list(),
                             AIT_anndata_path=NULL,
                             force=FALSE,
                             n_processors = 3,
                             normalization = "log2CPM",
                             tmp_dir = NULL,
                             user_precomp_stats_path=NULL,
                             user_query_markers_path=NULL){
  tryCatch(
    {
      # move to zzz try catch
      cell_type_mapper <- import("cell_type_mapper")
      temp_folder = tmp_dir

      if ((length(AIT_anndata$uns$hierarchical[[AIT_anndata$uns$mode]]) > 0) && force==FALSE) {
        stop(paste0(paste0("ERROR: mode provided '", AIT_anndata$uns$mode), 
        "' already exists, choose a new mode name or use force=TRUE to overwrite."))
      }

      if (is.null(temp_folder) || temp_folder == "") {
        temp_folder <- paste0("temp_folder_", format(Sys.time(), "%Y%m%d-%H%M%S"))
        temp_folder <- file.path(getwd(), temp_folder)
        dir.create(temp_folder)
      }

      # get an ordered list of taxonomy's hierarchy levels.
      taxonomy_hierarchy = get_hierarchy(AIT_anndata, hierarchy)

      # get file path to the AIT taxonomy (h5ad)
      taxonomy_anndata_path = get_anndata_path(AIT_anndata, AIT_anndata_path, temp_folder)

      # compute stats and save them to anndata
      precomp_stats_output_path = user_precomp_stats_path
      if(is.null(precomp_stats_output_path)) {
        precomp_stats_output_path = run_precomp_stats(taxonomy_anndata_path, n_processors, normalization, temp_folder, taxonomy_hierarchy)
      }
      AIT_anndata = save_precomp_stats_to_uns(taxonomy_anndata_path, precomp_stats_output_path, AIT_anndata$uns$mode)

      # compute query markers and save them to anndata
      query_markers_output_path = user_query_markers_path
      if(is.null(query_markers_output_path)) {
        ref_markers_file_path = run_reference_markers(precomp_stats_output_path, n_processors, temp_folder) 
        query_markers_output_path = run_query_markers(taxonomy_anndata_path, ref_markers_file_path, n_processors, temp_folder) 
      }
      AIT_anndata = save_query_markers_to_uns(AIT_anndata, query_markers_output_path) 
      
      AIT_anndata$write_h5ad(taxonomy_anndata_path)
    },
    error = function(e) {
      errorMessage <- conditionMessage(e)
      cat("Error message:", errorMessage, "\n")
    },
    finally = {
      # remove the temp folder is it was code generated
      if (is.null(tmp_dir) || tmp_dir == "") {
        print(paste("Deleting temp folder", temp_folder))
        unlink(temp_folder, recursive = TRUE)
      }
      else {
        # remove the files is they were code generated
        if(is.null(user_precomp_stats_path)) {
          file.remove(precomp_stats_output_path)
        }
        if(is.null(user_query_markers_path)) {
          file.remove(ref_markers_file_path)
          file.remove(query_markers_output_path)
        }
        if (is.null(AIT_anndata_path) && 
            is.null(AIT_anndata$uns$taxonomyDir) && 
            is.null(AIT_anndata$uns$title)) {
          file.remove(taxonomy_anndata_path)
        }
      }
      return(AIT_anndata)
    }
  )
}

#' Generates the precomputed_stats.h5 file, which is an HDF5 file that contains a serialization of the cell type taxonomy tree, 
#' as well as statistics about the reference cells that have been assigned to the cell types in the taxonomy. 
#'
#' @param anndata_path Local file path of the AIT reference taxonomy (h5ad file).
#' @param n_processors Number of independent worker processes to spin up.
#' @param normalization Normalization of the h5ad files; must be either 'raw' or 'log2CPM'.
#' @param temp_folder Temporary directory for writing out the hierarchical files (the code will clean these up after itself).
#' @param hierarchy List of term_set_labels in the reference taxonomy ordered from most gross to most fine.
#'
#' @return File path to the precompute stats file.
#'
#' @keywords internal
run_precomp_stats = function(anndata_path, n_processors, normalization, temp_folder, hierarchy) {

  ##
  temp_precomp_stats_name = paste0("precomp_stats_", format(Sys.time(), "%Y%m%d-%H%M%S"))
  precomp_stats_filename <- paste0(temp_precomp_stats_name, ".h5")
  precomp_stats_output_path <- file.path(temp_folder, precomp_stats_filename)

  precomp_stats_config <- list(
      'h5ad_path' = anndata_path,
      'n_processors' = n_processors,
      'normalization' = normalization,
      'tmp_dir' = temp_folder,
      'output_path' = precomp_stats_output_path,
      'hierarchy' = hierarchy
  )

  precomp_stats_runner <- cell_type_mapper$cli$precompute_stats_scrattch$PrecomputationScrattchRunner(
    args=c(),
    input_data=precomp_stats_config)
  precomp_stats_runner$run()

  return(precomp_stats_output_path)
}

#' Saves the contents of the precomputed_stats.h5 file to AIT -> uns -> hierarchical.
#'
#' @param anndata_path Local file path of the AIT reference taxonomy (h5ad file).
#' @param precomp_stats_output_path Local file path to the generated or user provided precomputed_stats.h5 file.
#'
#' @return AIT reference taxonomy with the precomupte stats saved in uns -> hierarchical.
#'
#' @keywords internal
save_precomp_stats_to_uns = function(anndata_path, precomp_stats_output_path, mode) {
  # save precomp stats to anndata h5ad file using cell_type_mapper's precomputed_stats_to_uns
  temp_precomp_stats_name = tools::file_path_sans_ext(basename(precomp_stats_output_path))
  cell_type_mapper$utils$output_utils$precomputed_stats_to_uns(
      h5ad_path=anndata_path,
      precomputed_stats_path=precomp_stats_output_path, 
      uns_key=temp_precomp_stats_name)

  ## load again b/c cell_type_mapper's precomputed_stats_to_uns saved precompstats to
  ## the .h5ad file, not AIT object.
  taxonomy_dir_path <- dirname(anndata_path)
  anndata_file_name <- basename(anndata_path)
  AIT_anndata = loadTaxonomy(taxonomyDir = taxonomy_dir_path, anndata_file=anndata_file_name)
  AIT_anndata$uns$mode = mode

  ## take the precomputed_stats from uns and save to uns$hierarchical$mode
  precomp_stats_json = AIT_anndata$uns[[temp_precomp_stats_name]]
  AIT_anndata$uns$hierarchical[[AIT_anndata$uns$mode]] <- list()
  AIT_anndata$uns$hierarchical[[AIT_anndata$uns$mode]][["precomp_stats"]] <- precomp_stats_json
  AIT_anndata$uns[[temp_precomp_stats_name]] <- NULL

  return(AIT_anndata)
}

#' Generates reference_marker.h5 file. The reference markers are, for every pair of leaf nodes in the cell type taxonomy tree, 
#' every gene that could conceivably be a marker gene for discriminating between those two cell types. 
#'
#' @param precomp_stats_output_path Local file path to the generated or user provided precomputed_stats.h5 file.
#' @param n_processors Number of independent worker processes to spin up.
#' @param temp_folder Temporary directory for writing out the hierarchical files (the code will clean these up after itself).
#'
#' @return File path of the reference marker file.
#'
#' @keywords internal
run_reference_markers = function(precomp_stats_output_path, n_processors, temp_folder) {
  ref_markers_config <- list(
      'n_processors' = n_processors,
      'precomputed_path_list' = list(precomp_stats_output_path),
      'output_dir' = temp_folder,
      'tmp_dir' = temp_folder
  )

  ref_markers_runner = cell_type_mapper$cli$reference_markers$ReferenceMarkerRunner(
      args=c(), 
      input_data=ref_markers_config)
  ref_markers_runner$run()

  ref_markers_file_path = file.path(temp_folder, "reference_markers.h5")
  return(ref_markers_file_path)
}

#' Generates query_markers.json file. This is done by converting reference_markers.h5 file into the final JSON lookup table of marker genes.
#'
#' @param anndata_path Local file path of the AIT reference taxonomy (h5ad file).
#' @param ref_markers_file_path Local file path to the generated reference_marker.h5 file.
#' @param n_processors Number of independent worker processes to spin up.
#' @param temp_folder Temporary directory for writing out the hierarchical files (the code will clean these up after itself).
#'
#' @return File path of the query markers file.
#'
#' @keywords internal
run_query_markers = function(anndata_path, ref_markers_file_path, n_processors, temp_folder) {
  query_markers_filename = paste0(paste0("query_markers_", format(Sys.time(), "%Y%m%d-%H%M%S")), ".json")
  query_markers_output_path = file.path(temp_folder, query_markers_filename)

  query_markers_config <- list(
      'query_path' = anndata_path,
      'reference_marker_path_list' = list(ref_markers_file_path),
      'n_processors' = n_processors,
      'output_path' = query_markers_output_path,
      'tmp_dir' = temp_folder
  )

  query_markers_runner = cell_type_mapper$cli$query_markers$QueryMarkerRunner(
      args=c(),
      input_data=query_markers_config)
  query_markers_runner$run()

  return(query_markers_output_path)
}

#' Saves the contents of the query_markers.json file to AIT -> uns -> hierarchical.
#'
#' @param AIT_anndata AIT reference taxonomy object where the file contents will be saved.
#' @param query_markers_output_path Local file path to the generated or user provided query_markers.json file.
#'
#' @return AIT reference taxonomy with the query markers saved in uns -> hierarchical.
#'
#' @keywords internal
save_query_markers_to_uns = function(AIT_anndata, query_markers_output_path) {
  ## extract query_markers data and serialize it
  query_markers_data = fromJSON(query_markers_output_path)
  serialized_query_markers = cell_type_mapper$utils$utils$clean_for_uns_serialization(query_markers_data)
  
  ## save serialized query_markers to uns$hierarchical$mode
  AIT_anndata$uns$hierarchical[[AIT_anndata$uns$mode]][["query_markers"]] <- serialized_query_markers

  return(AIT_anndata)
}

#' This function saves the AIT reference taxonomy to a temp folder as h5ad, if the provided file path is invalid.
#' @param AIT_anndata AIT reference taxonomy object.
#' @param anndata_path Local file path of the AIT reference taxonomy (h5ad file).
#' @param temp_folder Temporary directory for writing out temporary files (the code will clean these up after itself).
#' @return Local file path to the AIT reference taxonomy h5ad file.
#'
#' @keywords internal
get_anndata_path = function(AIT_anndata, anndata_path, temp_folder) {
  if (is.null(anndata_path) || !file.exists(anndata_path)){
    # Use AIT path stored in AIT.anndata$uns, if not null.
    if (!is.null(AIT_anndata$uns$taxonomyDir) && !is.null(AIT_anndata$uns$title)){
      anndata_path = file.path(AIT_anndata$uns$taxonomyDir, paste0(AIT_anndata$uns$title, ".h5ad"))
    }
    # Check if anndata path exists; if does not, write it out to temp - show WARNING.
    if (is.null(anndata_path) || !file.exists(anndata_path)) {
      print(paste0(paste("WARNING: INVALID FILE PATH, ERROR in AIT.anndata$uns taxonomyDir and taxonomyName:", anndata_path),
            ". Writing the AIT.anndata to temperary location, SAVE anndata or FIX path to OPTIMIZE this step."))
      anndata_filename <- paste0(paste0("temp_anndata_", format(Sys.time(), "%Y%m%d-%H%M%S")), ".h5ad")
      anndata_path <- file.path(temp_folder, anndata_filename)
      AIT_anndata$write_h5ad(anndata_path)
    }
  }
  return(anndata_path)
}

#' This function looks for a valid hierarchy list.
#' @param AIT_anndata AIT reference taxonomy object.
#' @param hierarchy List of term_set_labels in the reference taxonomy ordered from most gross to most fine.
#' @return Local file path to the AIT reference taxonomy h5ad file.
#'
#' @keywords internal
get_hierarchy = function(AIT_anndata, hierarchy) {
  if (length(hierarchy) == 0) {
    hierarchy = AIT_anndata$uns$hierarchy
  }
  else {
    AIT_anndata$uns$hierarchy = hierarchy 
  }
  if (is.null(hierarchy) || length(hierarchy) == 0) {
    stop(paste("Hierarchy does NOT exist in AIT.anndata$uns$hierarchy, please pass hierarchy list as a function parameter 'hierarchy=list()'."))
  }
  return(hierarchy)
}