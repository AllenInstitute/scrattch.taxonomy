#' @export
addMapMyCells = function(AIT_anndata,
                             hierarchy,
                             anndata_path=NULL,
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
        "' already exists, choose a new mode name or use force=TRUE to override."))
      }

      if(is.null(anndata_path)){
        anndata_path = file.path(AIT_anndata$uns$taxonomyDir, paste0(AIT_anndata$uns$title, ".h5ad"))
      }

      if (is.null(temp_folder) || temp_folder == "") {
        temp_folder <- paste0("temp_folder_", format(Sys.time(), "%Y%m%d-%H%M%S"))
        temp_folder <- file.path(getwd(), temp_folder)
        dir.create(temp_folder)
      }

      # compute stats and save them to anndata
      precomp_stats_output_path = user_precomp_stats_path
      if(is.null(precomp_stats_output_path)) {
        precomp_stats_output_path = run_precomp_stats(anndata_path, n_processors, normalization, temp_folder, hierarchy)
      }
      AIT_anndata = save_precopm_stats_to_uns(anndata_path, precomp_stats_output_path)

      # compute query markers and save them to anndata
      query_markers_output_path = user_query_markers_path
      if(is.null(query_markers_output_path)) {
        ref_markers_file_path = run_reference_markers(precomp_stats_output_path, n_processors, temp_folder) 
        query_markers_output_path = run_query_markers(anndata_path, ref_markers_file_path, n_processors, temp_folder) 
      }
      AIT_anndata = save_query_markers_to_uns(AIT_anndata, query_markers_output_path) 
      
      AIT_anndata$write_h5ad(anndata_path)
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
      }
      return(AIT_anndata)
    }
  )
}

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

#' @keywords internal
save_precopm_stats_to_uns = function(anndata_path, precomp_stats_output_path) {
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

  ## take the precomputed_stats from uns and save to uns$hierarchical$mode
  precomp_stats_json = AIT_anndata$uns[[temp_precomp_stats_name]]
  AIT_anndata$uns$hierarchical[[AIT_anndata$uns$mode]] <- list()
  AIT_anndata$uns$hierarchical[[AIT_anndata$uns$mode]][["precomp_stats"]] <- precomp_stats_json
  AIT_anndata$uns[[temp_precomp_stats_name]] <- NULL

  return(AIT_anndata)
}

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

#' @keywords internal
save_query_markers_to_uns = function(AIT_anndata, query_markers_output_path) {
  ## extract query_markers data and serialize it
  query_markers_data = fromJSON(query_markers_output_path)
  serialized_query_markers = cell_type_mapper$utils$utils$clean_for_uns_serialization(query_markers_data)
  
  ## save serialized query_markers to uns$hierarchical$mode
  AIT_anndata$uns$hierarchical[[AIT_anndata$uns$mode]][["query_markers"]] <- serialized_query_markers

  return(AIT_anndata)
}