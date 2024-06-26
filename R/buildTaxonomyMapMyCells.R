#' @export

buildHANNMapMyCells = function(AIT_anndata,
                             hierarchy,
                             anndata_path,
                             n_processors = 3,
                             normalization = "log2CPM",
                             tmp_dir = NULL){
  delete_tmp_folder <- FALSE
  tryCatch(
    {
      cell_type_mapper <- import("cell_type_mapper")
  
      if(missing(anndata_path)){
        anndata_path = file.path(AIT_anndata$uns$taxonomyDir, paste0(AIT_anndata$uns$taxonomyName, ".h5ad"))
      }

      if (missing(tmp_dir) || is.null(tmp_dir) || tmp_dir == "") {
        tmp_dir <- paste0("temp_folder_", format(Sys.time(), "%Y%m%d-%H%M%S"))
        tmp_dir <- file.path(getwd(), tmp_dir)
        dir.create(tmp_dir)
        delete_tmp_folder <- TRUE
      }

      # ========PRECOMPUTE STATS========
      precomp_stats_filename <- paste0(paste0("precomp_stats_", format(Sys.time(), "%Y%m%d-%H%M%S")), ".h5")
      precomp_stats_output_path <- file.path(tmp_dir, precomp_stats_filename)

      precomp_stats_config <- list(
          'h5ad_path' = anndata_path,
          'n_processors' = n_processors,
          'normalization' = normalization,
          'tmp_dir' = tmp_dir,
          'output_path' = precomp_stats_output_path,
          'hierarchy' = hierarchy
      )

      precomp_stats_runner <- cell_type_mapper$cli$precompute_stats_scrattch$PrecomputationScrattchRunner(
        args=c(),
        input_data=precomp_stats_config)
      precomp_stats_runner$run()

      # ADD TO UNS
      cell_type_mapper$utils$output_utils$precomputed_stats_to_uns(
          h5ad_path=anndata_path,
          precomputed_stats_path=precomp_stats_output_path, 
          uns_key="MapMyCells_HANN_precomp_stats")

      # ========REFERENCE MARKERS========
      ref_markers_config <- list(
          'n_processors' = n_processors,
          'precomputed_path_list' = list(precomp_stats_output_path),
          'output_dir' = tmp_dir,
          'tmp_dir' = tmp_dir
      )

      ref_markers_runner = cell_type_mapper$cli$reference_markers$ReferenceMarkerRunner(
          args=c(), 
          input_data=ref_markers_config)
      ref_markers_runner$run()


      # ========QUERY MARKERS========
      ref_markers_file_path = file.path(tmp_dir, "reference_markers.h5")
      query_markers_filename = paste0(paste0("query_markers_", format(Sys.time(), "%Y%m%d-%H%M%S")), ".json")
      query_markers_output_path = file.path(tmp_dir, query_markers_filename)

      query_markers_config <- list(
          'query_path' = anndata_path,
          'reference_marker_path_list' = list(ref_markers_file_path),
          'n_processors' = n_processors,
          'output_path' = query_markers_output_path,
          'tmp_dir' = tmp_dir
      )

      query_markers_runner = cell_type_mapper$cli$query_markers$QueryMarkerRunner(
          args=c(),
          input_data=query_markers_config)
      query_markers_runner$run()

      query_markers_data = fromJSON(query_markers_output_path)
      serialized_query_markers = cell_type_mapper$utils$utils$clean_for_uns_serialization(query_markers_data)
      query_markers_uns_key = "MapMyCells_HANN_query_markers"
      query_markers_uns = list()
      query_markers_uns[[query_markers_uns_key]] = serialized_query_markers

      cell_type_mapper$utils$anndata_utils$update_uns(
          h5ad_path=anndata_path,
          new_uns=query_markers_uns,
          clobber=FALSE
      )

      AIT_anndata = read_h5ad(anndata_path)
      return(AIT_anndata)
    },
    error = function(e) {
      errorMessage <- conditionMessage(e)
      cat("Error message:", errorMessage, "\n")
    },
    finally = {
      if (delete_tmp_folder) {
        print(paste("Deleting temp folder", tmp_dir))
        unlink(tmp_dir, recursive = TRUE)
      }
    }
  )
}