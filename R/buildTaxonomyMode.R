#' Create a taxonomy mode for mapping to taxonomy subset (or child)
#'
#' This function creates a new mode for mapping to a subset of a taxonomy. This includes filtering cells, setting new variable genes and/or marker genes, calculating mapping statistics, 
#'
#' @param AIT.anndata A reference taxonomy anndata object.
#' @param mode.name A name to identify the new taxonomy version.
#' @param retain.cells A boolean vector of length [number of cells] indicating which cells should be retained (TRUE) or filtered (FALSE) -OR- a character vector with sample names indicating which cells should be retained. Default is to retain the cells included in "stadard" mode.
#' @param retain.clusters A character vector with cluster names (e.g., values in the "cluster_id" column) indicating which clusters should be retained. Default is to retain all clusters with at least 2 retained cells (clusters with exactly 1 cell can cause some functions to crash).
#' @param subsample The number of cells to retain per cluster (default = 100). Note that subsampling happens AFTER retain.cells and retail.clusters filtering.
#' @param highly_variable_genes Set of features defined as highly variable genes. Provide either as a named list of vectors, or as a single vector (in which case the name "highly_variable_genes_[mode.name]" will be used). Optional input, but for proper mapping we recommend including either highly_variable_genes or marker_genes. If neither are provided, standard mode markers will be used for mapping algorithms based on these gene lists and may cause problems.
#' @param marker_genes Set of features defined as marker genes. Provide either as a named list of vectors, or as a single vector (in which case the name "marker_genes_[mode.name]" will be used).
#' @param embeddings Dimensionality reduction coordinate data.frame with 2 columns. Rownames must be equal to colnames of counts.  Either provide as a named list or as a single data.frame (in which case the name "default_[mode.name]" will be used). Optional - if not provided the relevant subset of the default standard embedding will be used.
#' @param add.dendrogram.markers If TRUE (default=FALSE), will also add dendrogram markers to prep the taxonomy for tree mapping
#' @param addMapMyCells If TRUE (default), will also prep this mode of the taxonomy for hierarchical MapMyCells mapping
#' @param overwrite If [mode.name] already exists, should it be overwritten (default = FALSE)
#' @param ... Additional variables to be passed to `addDendrogramMarkers`
#' 
#' @import scrattch.hicat
#' @import reticulate
#'
#' @return AIT.anndata An updated AIT.anndata variable with new content addded for the relevant mode.name (also called "Child taxonomy").
#'
#' @export
buildTaxonomyMode = function(AIT.anndata,
                             mode.name,
                             retain.cells = NULL,
                             retain.clusters = NULL,
                             subsample = 100,
                             highly_variable_genes = NULL,
                             marker_genes = NULL,
                             embeddings = NULL,
                             add.dendrogram.markers = FALSE,
                             addMapMyCells = TRUE,
                             overwrite = FALSE,
                             ...){
  
  ## FIRST DO SOME CHECKS OF THE INPUT VARIABLES
  
  ## Ensure filtering mode doesn't already exist
  if(!is.character(mode.name)) stop("mode.name must be a character string.")
  mode.name <- mode.name[1]
  if(mode.name %in% names(AIT.anndata$uns$filter)){ 
    if(!overwrite) stop(print(paste0("Mode ", mode.name, " already exists in the taxonomy. To rebuild this mode, set 'overwrite=TRUE'.")))
    paste0("Mode ", mode.name, " already exists in the taxonomy. You will be overwriting the previous mode files.")
  }
  
  ## If highly_variable_genes and/or marker genes are provided as vectors (as with previous versions of scrattch.taxonomy, convert to named lists)
  if(is.character(highly_variable_genes)){
    highly_variable_genes <- list(highly_variable_genes = highly_variable_genes)
    names(highly_variable_genes) <- paste0("highly_variable_genes_",mode.name)
  }
  if(is.character(marker_genes)){
    marker_genes <- list(marker_genes = marker_genes)
    names(marker_genes) <- paste0("marker_genes_",mode.name)
  }
  
  ## Determine the cluster column and warn if not "cluster_id"
  hierarchy = AIT.anndata$uns$hierarchy
  hierarchy = hierarchy[order(unlist(hierarchy))]
  if(is.null(hierarchy)) stop("Hierarchy must be included in the standard AIT mode in proper format to create a mode.  Please run checkTaxonomy().")
  celltypeColumn = names(hierarchy)[length(hierarchy)][[1]]
  if(celltypeColumn!="cluster_id") warning("AIT schema requires clusters to be in 'cluster_id' slot. We recommend calling the finest level of the hierarch as 'cluster_id'.")
  
  
  ## NOW FILTER THE CELLS AND CLUSTERS, AND SUBSAMPLE
  
  sample.vector  = row.names(AIT.anndata$obs)
  cluster.vector = as.character(AIT.anndata$obs[,celltypeColumn])
  retain = rep(TRUE, length(cluster.vector))
  if(!is.null(retain.cells)){
    if(is.logical(retain.cells)){
      retain = retain.cells
    } else{
      retain = sample.vector %in% retain.cells
    }
  }
  if(!is.null(retain.clusters)){
    if(is.logical(retain.clusters)){
      stop("If provided, cluster vector must be a character vector of cluster names.")
    } else{
      retain = retain & (cluster.vector %in% retain.clusters)
    }
  }
  if(sum(retain)<2) stop("No cells are retained! Double check filtering options in input.")
  retain.subsample = subsampleCells(cluster.vector[retain],subsample)
  retain[retain] = retain.subsample
  mode.clusters  = sort(unique(cluster.vector[retain]))
  if(length(mode.clusters)<2) stop("At least two clusters are needed for a proper taxonomy. Double check filtering options in input.")
  AIT.anndata$uns$filter[[mode.name]] = !retain
  AIT.anndata$uns$stats[[mode.name]][["clusters"]] = mode.clusters

  
  ## Filter stats files if they exist
  #AIT.anndata$varm[[paste0("cluster_id_median_expr_",mode.name)]] =  AIT.anndata$varm[["cluster_id_median_expr_standard"]][,mode.clusters]
  #  THIS IS REDUNDANT TO SAVE, SO WE WON'T.  IT CAN BE SUBSET AS NEEDED
  
  
  ## MODIFY THE DENDROGRAM AND SAVE IN THE ANNDATA
  
  if(!is.null(AIT.anndata$uns$dend)){
    dend = json_to_dend(AIT.anndata$uns$dend[["standard"]])
    dend = prune(dend, setdiff(labels(dend), mode.clusters))
    AIT.anndata$uns$dend[[mode.name]] = toJSON(dend_to_json(dend))
  }

  
  ## Update markers after pruning
  if ((add.dendrogram.markers)&(!is.null(AIT.anndata$uns$dend))){
    AIT.anndata = addDendrogramMarkers(AIT.anndata, mode=mode.name, taxonomyDir=taxonomyDir, ...)
    # The reference probability matrix for the subsetted taxonomy is defined and outputted in this function as well
    # $memb[[mode.name]]
    # ...$memb.ref,
    # ...$map.df.ref
  }
  
  
  ## Add MapMyCells (hierarchical mapping) functionality, if requested
  if(addMapMyCells) {
    current.mode = AIT.anndata$uns$mode
    AIT.anndata  = mappingMode(AIT.anndata, mode=mode.name)
    AIT.anndata  = addMapMyCells(AIT.anndata, hierarchy, force=TRUE)
    AIT.anndata  = mappingMode(AIT.anndata, mode=current.mode)
    # NOTE: statistics used to be saved in 'uns', but are being moved to 'varm'
  }
  
  ## If not yet written in add.dendrogram.markers or addMapMyCells, write h5ad
  if(!(addMapMyCells|add.dendrogram.markers)){
    AIT.anndata$write_h5ad(file.path(AIT.anndata$uns$taxonomyDir, paste0(AIT.anndata$uns$title, ".h5ad")))
  }
  
  ##
  return(AIT.anndata)

}
