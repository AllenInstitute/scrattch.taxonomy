#' This function builds the minimum files required for Shiny
#'
#' @param counts A count matrix in sparse format: dgCMatrix.
#' @param meta.data Meta.data corresponding to count matrix. Rownames must be equal to colnames of counts. "clusters" must be provided (see celltypeColumn and notes).
#' @param feature.set Set of feature used to calculate dendrogram. Typically highly variable and/or marker genes.
#' @param umap.coords Dimensionality reduction coordiant data.frame with 2 columns. Rownames must be equal to colnames of counts.
#' @param subsample The number of cells to retain per cluster (default = 2000)
#' @param taxonomyDir The location to save Shiny objects, e.g. "/allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/10x_seq/NHP_BG_20220104/"
#' @param taxonomyTitle The file name to assign for the Taxonomy h5ad.
#' @param celltypeColumn Column name corresponding to where the clusters are located (default="cluster")
#' @param hierarchy List of term_set_labels in the Taxonomy ordered from most gross to most fine (e.g., subclass, supertype, neighborhood, class).
#' @param cluster_colors An optional named character vector where the values correspond to colors and the names correspond to celltypes in celltypeColumn.  If this vector is incomplete, a warning is thrown and it is ignored. cluster_colors can also be provided in the metadata (see notes)
#' @param cluster_stats A matrix of median gene expression by cluster. Cluster names must exactly match meta.data$cluster.
#' @param dend Existing dendrogram associated with this taxonomy (e.g., one calculated elsewhere).  If NULL (default) a new dendrogram will be calculated based on the input `feature.set`
#' @param reorder.dendrogram Should dendogram attempt to match a preset order? (Default = FALSE).  If TRUE, the dendrogram attempts to match the celltype factor order as closely as possible (if celltype is a character vector rather than a factor, this will sort clusters alphabetically, which is not ideal).
#' @param add.dendrogram.markers If TRUE (default), will also add dendrogram markers to prep the taxonomy for tree mapping
#' @param addMapMyCells If TRUE (default), will also prep this taxonomy for hierarchical mapping
#' @param ... Additional variables to be passed to `addDendrogramMarkers`
#' 
#' Notes: Precomputed clusters must be provided.  In the anndata object these will be stored using the term "cluster".  If celltypeColumn is anything other than cluster, then any existing "cluster" column will be overwritten by celltypeColumn.  Values can be provided without colors and ids (e.g., "cluster") or with them (e.g., "cluster_label" + "cluster_color" + "cluster_id").  In this case cluster_colors is ignored and colors are taken directly from the metadata.  Cluster_id's will be overwritten to match dendrogram order.
#' 
#' @import scrattch.hicat 
#' @import feather
#' @import tibble
#' @import dplyr
#' @import Matrix
#' @import pvclust
#' @import anndata
#' @import reticulate
#'
#' @return AIT anndata object in the specified format (only if return.anndata=TRUE)
#'
#' @export
buildTaxonomy = function(counts,
                         meta.data,
                         feature.set,
                         umap.coords,
                         tpm = NULL,
                         subsample=2000,
                         taxonomyDir = getwd(),
                         taxonomyTitle = "AI_taxonomy",
                         celltypeColumn = "cluster",
                         hierarchy = list(),
                         cluster_colors = NULL,
                         cluster_stats = NULL,
                         dend = NULL,
                         reorder.dendrogram = FALSE,
                         add.dendrogram.markers = TRUE,
                         addMapMyCells = TRUE,
                         ...){

  ## Sanity check and cleaning of parameters
  clean.params = .checkBuildTaxonomyParams(counts, 
                                            meta.data, 
                                            feature.set, 
                                            umap.coords, 
                                            taxonomyDir, 
                                            taxonomyTitle, 
                                            celltypeColumn, 
                                            cluster_colors, 
                                            cluster_stats,
                                            dend)
  meta.data = clean.params$meta.data; celltypeColumn = clean.params$celltypeColumn

  ## ----------
  ## Run auto_annotate, this changes cell_id to cell_id.
  meta.data$cell_id = colnames(counts)
  meta.data = .formatMetadata(meta.data, cluster_colors)
  
  ## ----------
  ## Subsample nuclei per cluster, max of subsample cells per cluster
  kpClusters <- rep(TRUE,length(meta.data$cluster_label))
  
  if(!is.null(dend))
    kpClusters <- is.element(meta.data$cluster_label, labels(dend)) # exclude missing clusters, if any

  if((subsample > 0) & (subsample < Inf)){
      print("===== Subsampling cells =====")
      kpSub = colnames(counts)[subsampleCells(meta.data$cluster_label, subsample) & kpClusters]
  }else{
      print("===== No subsampling of cells =====")
      kpSub = colnames(counts)[kpClusters]
  }

  ## ----------
  ## Ensure meta.data is labeled with cell identifiers
  rownames(meta.data) = meta.data$cell_id

  ## Gather clusters
  clustersUse = unique(meta.data$cluster_label)

  ## Define clusterInfo, which is used to convert cell types to subclass / neighborhood / class
  clusterInfo = as.data.frame(meta.data[match(clustersUse, meta.data$cluster_label),])

  ## Computing TPM matrix
  if(is.null(tpm)){
    print("===== Normalizing count matrix to log2(CPM) =====")
    tpm = as(scrattch.bigcat::logCPM(counts), "dgCMatrix")
  }else{
    print("===== User provided TPM matrix =====")
  }

  ## Compute cluster medians if needed
  if(is.null(cluster_stats)){
    print("===== Computing median expr. at taxonomy leaves =====")
    ## Get cluster medians
    cluster   = meta.data$cluster_label; names(cluster) = rownames(meta.data)
    medianExpr = scrattch.bigcat::get_cl_medians(tpm, cluster) 
  }else{
    print("===== User provided median expr. at taxonomy leaves =====")
    medianExpr = cluster_stats
    ## Sanity check user input
  }

  ##  Build the dendrogram
  print("===== Building dendrogram =====")
  if(!is.null(dend)){
    print("...using provided dendrogram.")
    # FOR FUTURE UPDATE: should check here whether dendrogram colors match what is in meta-data.
  } else {
    ## Define the cluster info 
    unique.meta.data = meta.data %>% distinct(cluster_id, 
                                              cluster_label, 
                                              cluster_color)
    rownames(unique.meta.data) = unique.meta.data$cluster_label
    ## Dendrogram parameters and gene sets
    use.color = setNames(unique.meta.data$cluster_color, unique.meta.data$cluster_label)[colnames(medianExpr)]
    if(reorder.dendrogram){
      l.rank = setNames(meta.data$cluster_id[match(unique.meta.data$cluster_label, meta.data$cluster_label)], unique.meta.data$cluster_label)
      l.rank = sort(l.rank)
    }else{ l.rank    = NULL }
    ## Build the dendrogram
    invisible(capture.output({  # Avoid printing lots of numbers to the screen
      dend.result = build_dend(
        cl.dat  = medianExpr[feature.set,],
        cl.cor  = NULL,
        l.color = use.color,
        l.rank  = l.rank, 
        nboot   = 1,
        ncores  = 1)
    }))
    dend = dend.result$dend
    print("...dendrogram built.")
  }

  ##
  print("===== Building taxonomy anndata =====")
  AIT.anndata = AnnData(
    X = Matrix::t(tpm), ## logCPM ensured to be in sparse column format
    raw = list(
      X = as(Matrix::t(counts), "dgCMatrix"),
      var = data.frame("gene" = rownames(counts))
    ), ## Store counts matrix
    obs = meta.data,
    var = data.frame("gene" = rownames(counts), 
                     "highly_variable_genes" = rownames(counts) %in% feature.set,          # retain original name for back-compatibility
                     "highly_variable_genes_standard" = rownames(counts) %in% feature.set, # copy to "_standard" for standard mode
                     row.names=rownames(counts)),
    obsm = list(
      X_umap = umap.coords ## A data frame with cell_id, and 2D coordinates for umap (or comparable) representation(s)
    ),
    uns = list(
      dend        = list("standard" = toJSON(dend_to_json(dend))), ## JSON dendrogram
      filter      = list("standard" = !(colnames(counts) %in% kpSub)), ## Filtered cells
      QC_markers  = list("standard" = list()), ## Standard will hold de.genes for dendrogram, we should rename this uns field.
      stats   = list("standard" = list("medianExpr" = medianExpr, "features" = rownames(medianExpr), "clusters" = colnames(medianExpr))), ## Median expression by cluster
      hierarchical = list("standard" = list()),
      mode = "standard", ## Default mode to standard
      cellSet = colnames(counts),
      clustersUse = clustersUse,
      clusterInfo = clusterInfo,
      title = taxonomyTitle,
      hierarchy = hierarchy,
      taxonomyDir = file.path(normalizePath(taxonomyDir), leading_string="/") ## Normalize path in case where user doesn't provide absolute path.
    )
  )
  
  ## Add dendrogram markers and membership tables, if requested
  if(add.dendrogram.markers){
    print("===== Adding dendrogram markers and membership tables for tree mapping and patch-seq QC =====")
    AIT.anndata = addDendrogramMarkers(AIT.anndata, mode="standard", taxonomyDir=taxonomyDir, ...)
  # The reference probability matrix for the subsetted taxonomy is defined and outputted in this function as well
  # $memb[[mode.name]]
  # ...$memb.ref,
  # ...$map.df.ref
  }
  
  ## Add MapMyCells (hierarchical mapping) functionality, if requested
  if(addMapMyCells) {
    if(is.character(hierarchy)) hierarchy <- as.list(hierarchy)
    if((length(hierarchy)==0)|(sum(class(hierarchy)=="list")<1)){
      warning("hierarchy must be a list of term_set_labels in the reference taxonomy ordered from most gross to most fine included in AIT_anndata or provided separately. Since this is NOT the case, addMapMyCells is being skipped")
    } else{
      print("===== Adding MapMyCells (hierarchical mapping) functionality =====")
      AIT.anndata = mappingMode(AIT.anndata, mode="standard")
      AIT.anndata = addMapMyCells(AIT.anndata, hierarchy, force=TRUE)
    }
  }

  print("===== Writing taxonomy anndata =====")
  AIT.anndata$write_h5ad(file.path(AIT.anndata$uns$taxonomyDir, paste0(AIT.anndata$uns$title, ".h5ad")))
  
  ## Return the anndata object
  return(AIT.anndata)
}