#' This function builds the minimum files required for Shiny
#' Precomputed clusters must be provided.  In the anndata object these will be stored using the term "cluster".  If hierarchy[[-1]] is anything other than cluster, then any existing "cluster" column will be overwritten by hierarchy[[-1]].  Values can be provided without colors and ids (e.g., "cluster") or with them (e.g., "cluster_label" + "cluster_color" + "cluster_id").  In this case cluster_colors is ignored and colors are taken directly from the metadata.  Cluster_id's will be overwritten to match dendrogram order.  (NOTE: Some functionality with metadata colors is still under development.)
#'
#' @param meta.data Meta.data corresponding to count matrix. Rownames must be equal to colnames of counts. "clusters" must be provided (see hierarchy[[-1]] and notes).
#' @param title The file name to assign for the Taxonomy h5ad (default="AIT"; recommended to create your own title!).
#' @param counts A count matrix in sparse format: dgCMatrix. buildTaxonomy can work with count matrices that have cells are rows or columns, so long as counts has both row names AND column names.
#' @param highly_variable_genes Set of features defined as highly variable genes OR a number of binary genes to calculate (we recommend ~1000 - ~5000, for <100 to ~5000 cell types). If a feature list is provided, provide either as a named list of vectors, or as a single vector (in which case the name "highly_variable_genes_standard" will be used). "highly_variable_genes_standard" will also be used for calculated variable genes. Optional input, but for proper mapping we strongly recommend including either highly_variable_genes or marker_genes. 
#' @param marker_genes Set of features defined as marker genes. Provide either as a named list of vectors, or as a single vector (in which case the name "marker_genes_[mode.name]" will be used).
#' @param ensembl_id A vector of ensembl ids corresponding to the gene symbols in counts.
#' @param gene.meta.data Either NULL (default) or a data frame of additional gene information to include in the var component of anndata
#' @param cluster_stats A matrix of median gene expression by cluster. Cluster names must exactly match meta.data$cluster.  If provided, will get saved to "varm$cluster_id_median_expr_[mode]"
#' @param embeddings Dimensionality reduction coordinate data.frame with 2 columns or a string with the column name for marker_genes or variable_genes from which a UMAP should be calculated. If coordinates are provided, rownames must be equal to colnames of counts.  Either provide as a named list or as a single data.frame (in which case the name "default_standard" will be used). embeddings are not required, but inclusion of at least one embedding is strongly recommended.#' 
#' @param number.of.pcs Number of principle components to use for calculating UMAP coordinates (default=30). This is only used in embeddings corresponds to a variable gene column from which a UMAP should be calculated.
#' @param dend Existing dendrogram associated with this taxonomy (e.g., one calculated elsewhere). If provided, dend must be a BINARY tree. Can also input a string with the column name for marker_genes or variable_genes from which a dendrogram should be calculated. If NULL or if the function can't figure out what gene set you want, no dendrogram will be calculated. Failure to define a dendrogram will prevent some mapping algorithms from working properly! The default is to build a dendrogram using the first highly_variable_genes or marker_genes set.
#' @param taxonomyDir The location to save Shiny objects, e.g. "/allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/10x_seq/NHP_BG_20220104/"
#' @param hierarchy List of term_set_labels in the Taxonomy ordered from most gross to most fine (e.g., neighborhood, class, subclass, supertype).
#' @param subsample The number of cells to retain per cluster (default = 2000)
#' @param cluster_colors An optional named character vector where the values correspond to colors and the names correspond to celltypes in hierarchy[[-1]].  If this vector is incomplete, a warning is thrown and it is ignored. cluster_colors can also be provided in the metadata (see notes)
#' @param default_embedding A string indicating which embedding to use for calculations.  Default (NULL) is to take the first one provided in embeddings.
#' @param uns.variables If provided, a list of additional variables to be included in the uns.  See Notes for schema variables not otherwise accounted for.
#' @param reorder.dendrogram Should dendogram attempt to match a preset order? (Default = FALSE).  If TRUE, the dendrogram attempts to match the celltype factor order as closely as possible (if celltype is a character vector rather than a factor, this will sort clusters alphabetically, which is not ideal).
#' @param add.dendrogram.markers If TRUE (default=FALSE), will also add dendrogram markers to prep the taxonomy for tree mapping
#' @param addMapMyCells If TRUE (default), will also prep this taxonomy for hierarchical mapping
#' @param check.taxonomy Should the taxonomy be checked to see if it follows the AIT schema (default=TRUE)
#' @param print.messages If check.taxonomy occurs, should any messages be written to the screen in addition to the log file (default=TRUE)
#' @param ... Additional variables to be passed to `addDendrogramMarkers`
#' 
#' *Additional uns.variables*:
#' * dataset_purl: Link to molecular data if not present in X or raw.X.
#' * batch_condition: Keys defining batches for normalization/integration algorithms. Used for cellxgene.
#' * reference_genome: Reference genome used to align molecular measurements.
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
buildTaxonomy = function(title="AIT",
                         meta.data,
                         hierarchy,
                         ## X
                         counts = NULL, ## matrix, sparse.
                         normalized.expr = NULL, ## matrix, sparse.
                         ## var
                         highly_variable_genes = NULL, ## named list
                         marker_genes = NULL, ## named list
                         ensembl_id = NULL,
                         gene.meta.data = NULL,
                         ## varm
                         cluster_stats = NULL,
                         ## obsm
                         embeddings = NULL,
                         number.of.pcs = 30, 
                         ## uns
                         dend = NULL,  # This default is to choose a sensible gene set
                         taxonomyDir = getwd(),
                         cluster_colors = NULL, ## @Jeremy please handle this.
                         default_embedding = NULL,
                         uns.variables = list(),  ## @Nelson: this is my inelegant way of handling Additional uns.variables (see notes above... feel free to edit)
                         ## Additional parameters
                         subsample=2000,
                         reorder.dendrogram = FALSE,
                         add.dendrogram.markers = FALSE,
                         addMapMyCells = TRUE,
                         check.taxonomy = TRUE,
                         print.messages = TRUE,
                         ...){
  
  ## Setting NA values to 0 the meta.data; otherwise write.h5ad will crash
  if(sum(is.na(meta.data))>0) {
    warning("Attempting to set 'NA' meta.data values to 0; otherwise, write.h5ad will crash. If error occurs at the'===== Writing taxonomy anndata =====', this step failed and NA values must be removed prior to running buildTaxonomy.")
    meta.data[is.na(meta.data)] = 0 
  }
  
  ## Prime features.dendrogram (this may not be needed)
  features.dendrogram = NULL

  ## Ensure that hierarchy is a named list with ascending order to hierarchy, e.g. ["Class"=0, "Subclass"=1, "cluster"=2]
  ##   --- Also apply some checks for previous versions
  if(is.character(hierarchy))
    hierarchy = as.list(hierarchy)
  if(class(hierarchy[[1]])=="character"){
    ordered_hierarchy = setNames(as.list(seq_along(hierarchy)-1), unlist(hierarchy))
  } else {
    ordered_hierarchy = hierarchy
  }
  if(!all(names(ordered_hierarchy) %in% colnames(meta.data))) {
    stop("Hierarchy values must all be included as column names in the metadata.")
  }
  hierarchy = ordered_hierarchy
  
  ## Pull the finest level cell type column
  celltypeColumn = names(hierarchy)[length(hierarchy)][[1]]
  if(celltypeColumn!="cluster_id") warning("AIT schema requires clusters to be in 'cluster_id' slot. We recommend calling the finest level of the hierarchy as 'cluster_id'.")

  ## Transpose counts and convert to dgCMatrix if needed
  if(dim(counts)[2]==dim(meta.data)[1]){
    t.counts <- as(counts, "dgCMatrix")
    counts   <- Matrix::t(counts)
  }
  if((!is.null(counts))&(!("dgCMatrix" %in% as.character(class(taxonomy.counts)))))
    counts <- as(counts, "dgCMatrix")
  
  ## Deal with character dendrograms
  if(is.character(dend)){
    features.dendrogram = dend
    dend = NULL
  }
  
  ## Sanity check and cleaning of parameters
  clean.params = .checkBuildTaxonomyParams(counts, 
                                           normalized.expr,
                                           meta.data, 
                                           highly_variable_genes,
                                           marker_genes,
                                           embeddings, 
                                           celltypeColumn,
                                           cluster_stats,
                                           taxonomyDir, 
                                           title, 
                                           dend)

  ## ----------
  ## Subsample nuclei per cluster, max of subsample cells per cluster
  kpSub = subsample_taxonomy(meta.data[[celltypeColumn]], rownames(meta.data), dend, subsample)

  ## Define clusterInfo, which is used to convert cell types to subclass / neighborhood / class
  clusterInfo = meta.data %>% distinct(.data[[celltypeColumn]], .keep_all = TRUE)
  rownames(clusterInfo) = clusterInfo[[celltypeColumn]]

  ## Computing TPM matrix if count matrix exists
  if(!is.null(counts) & is.null(normalized.expr)){
    print("===== Normalizing count matrix to log2(CPM+1) =====")
    normalized.expr = log2CPM_byRow(counts)
  }else{
    print("===== No provided count matrix. Skipping TPM calculation. =====")
  }

  ## Compute cluster medians if needed
  if(!is.null(normalized.expr)){
    if(is.null(cluster_stats)){
      print("===== Computing median expr. at taxonomy leaves =====")
      ## Get cluster medians
      cluster   = meta.data[[celltypeColumn]]; names(cluster) = rownames(meta.data)
      t.normalized.expr <- Matrix::t(normalized.expr)  # We still need to transpose this matrix since get_cl_medians is SUPER complicated code
      t.normalized.expr <- as(t.normalized.expr,"dgCMatrix")
      cluster_stats = scrattch.bigcat::get_cl_medians(t.normalized.expr, cluster)
    }
  }
  
  ## Calculate binary genes as highly variable genes, if requested
  # If numeric, calculating the top N binary genes and save as a list first
  if(is.numeric(highly_variable_genes)){
    highly_variable_genes <- round(max(min(dim(counts)[2],highly_variable_genes[1]),100))  # Make sure the numeric value is legal
    print(paste0("===== Computing ",highly_variable_genes," binary genes as highly variable genes. ====="))
    if(!exists("t.counts")){
      t.counts <- Matrix::t(counts)         # Still need to take the transpose here because get_cl_props is really complicated code
      t.counts <- as(t.counts,"dgCMatrix")
    }
    binary.genes          <- top_binary_genes(t.counts, meta.data[[celltypeColumn]], highly_variable_genes)
    highly_variable_genes <- list("highly_variable_genes_standard" = binary.genes)
  }
  
  ## If highly_variable_genes and/or marker genes are provided as vectors (as with previous versions of scrattch.taxonomy, convert to named lists)
  if(is.character(highly_variable_genes)){
    highly_variable_genes <- list(highly_variable_genes_standard = highly_variable_genes)
  }
  if(is.character(marker_genes)){
    marker_genes <- list(marker_genes_standard = marker_genes)
  }
  
  ## Build the AIT object
  print("===== Building taxonomy anndata =====")
  AIT.anndata = AnnData(
    X = if(!is.null(normalized.expr)) (normalized.expr) else NULL, ## logCPM will already be in sparse column format, if provided
    raw = if(!is.null(counts)) list(X = counts, var = data.frame("gene" = colnames(counts))) else NULL, ## Store counts if provided
    obs = meta.data,
    var = if(!is.null(counts))
                    data.frame("gene" = colnames(counts), 
                               row.names = colnames(counts))
          else data.frame(),
    varm = list(
      "cluster_id_median_expr_standard" = cluster_stats ## Median expression by cluster
    ),
    obsm = list(  ## A data frame with cell_id, and 2D coordinates for umap (or comparable) representation(s)
    ),  
    uns = list(
      filter     = list("standard" = !(rownames(counts) %in% kpSub)), ## Filtered cells
      mapmycells = list("standard" = list()),
      mode = "standard", ## Default mode to standard
      cellSet = rownames(meta.data),
      cluster_info = clusterInfo,
      clusterStatsColumns = list("standard" = colnames(cluster_stats)),
      title = title,
      hierarchy = hierarchy,
      taxonomyDir = file.path(normalizePath(taxonomyDir), leading_string="/"), ## Normalize path in case where user doesn't provide absolute path.
      schema_version = as.character(packageVersion("scrattch.taxonomy"))
    )
  )
  
  print(".... Adjusting a few variables")
  
  ## Ensure the hierarchy is correctly ordered and not alphabetical (this shouldn't be necessary, but is, and still might not take)
  AIT.anndata$uns$hierarchy <- AIT.anndata$uns$hierarchy[order(as.numeric(AIT.anndata$uns$hierarchy))]

  ## highly_variable_genes is a data.frame with gene names in rows and various sets in columns
  if(!is.null(highly_variable_genes)){
    for(feature_set in names(highly_variable_genes)){
      AIT.anndata$var[[feature_set]] = rownames(AIT.anndata$var) %in% highly_variable_genes[[feature_set]]
    }
  }

  ## marker_genes is a data.frame with gene names in rows and various sets in columns
  if(!is.null(marker_genes)){
    for(feature_set in names(marker_genes)){
      AIT.anndata$var[[feature_set]] = rownames(AIT.anndata$var) %in% marker_genes[[feature_set]]
    }
  }
  
  ## Determine feature set for the dendrogram and then define the dendrogram, if requested
  if((is.null(dend))&(!((all(is.na(dend))&(length(dend[1])==1))))){
    print(".... determining gene set for calculating dendrogram.")
    dend = features.dendrogram
    features.dendrogram = NULL
    if(is.character(dend)){
      if(length(dend)>1){
        features.dendrogram <- intersect(features.dendrogram,colnames(AIT.anndata$raw$X))
        print(".... Setting dendrogram features as the set of input genes overlapping with available genes.")
      }
      if(length(dend)==1){
        dend <- intersect(dend,colnames(AIT.anndata$var))
        print(paste("Setting dendrogram features as",dend))
        if(!is.null(dend))
          features.dendrogram <- AIT.anndata$var[[dend]]
      }
    }
    if(length(features.dendrogram)<10) 
      features.dendrogram = NULL  # In this case, we cannot figure out which genes to use so all genes will be used
    dend=NULL
  }
  
  ##  Build the dendrogram
  print("===== Building dendrogram =====")
  if(!is.null(dend)){
    print("...using provided dendrogram.")
    # FOR FUTURE UPDATE: should check here whether dendrogram colors match what is in meta-data.
  } else {
    ## Dendrogram parameters and gene sets
    # use.color = setNames(clusterInfo$cluster_color, clusterInfo[[celltypeColumn]])[colnames(cluster_stats)] @Jeremy please handle this.
    if(reorder.dendrogram){
      l.rank = setNames(meta.data[[celltypeColumn]][match(clusterInfo[[celltypeColumn]], meta.data[[celltypeColumn]])], clusterInfo[[celltypeColumn]])
      l.rank = sort(l.rank)
    } else { 
      l.rank = NULL 
    }
    ## Figure out feature set to use 
    if(is.null(features.dendrogram)){
      feature.set = rownames(cluster_stats)
    } else {
      feature.set = features.dendrogram
    }
    ## Build the dendrogram
    if(all(is.na(dend))&(length(dend)==1)){
      dend = NULL
      print("...dendrogram IS NOT BUILT. This may cause issues for mapping!")
      warning("Dendrogram IS NOT BUILT. This may cause issues for mapping!")
    } else {
      invisible(capture.output({  # Avoid printing lots of numbers to the screen
        dend.result = build_dend(
          cl.dat  = cluster_stats[feature.set,],
          cl.cor  = NULL,
          # l.color = use.color, @Jeremy please handle this.
          l.rank  = l.rank, 
          nboot   = 1,
          ncores  = 1)
      }))
      dend <- dend.result$dend
      print("...dendrogram built.")
    }
  }
  
  ## Add the dendrogram to the AIT file, if available
  if(!is.null(dend))
    AIT.anndata$uns$dend = list("standard" = toJSON(dend_to_json(dend))) ## JSON dendrogram
  
  
  ## Add ensembl_id into AIT object
  if(!is.null(ensembl_id)){
    if(dim(AIT.anndata$var)[1]==length(ensembl_id)){
      print("===== Adding provided ensembl_id vector. =====")
      AIT.anndata$var$ensembl_id = as.character(ensembl_id)
    } else {
      print("===== Provided ensembl_id vector of different length from gene length in var. Skipping adding ensembl_ids. =====")
    }
  }
  
  ## Add other gene information fields, if provided
  if(!is.null(gene.meta.data)){
    if(dim(AIT.anndata$var)[1]==dim(gene.meta.data)[1]){
      print("===== Adding provided gene.meta.data. =====")
      AIT.anndata$var = cbind(AIT.anndata$var,gene.meta.data)
    } else {
      print("===== Provided gene.meta.data data frame has different row dimensions from other gene fields. Skipping. =====")
    }
  }
  
  ## Calculate embeddings if requested
  if(is.character(embeddings)){
    embeddings <- intersect(embeddings,colnames(AIT.anndata$var))
    if((!is.null(embeddings))&(!is.null(normalized.expr))){
      embeddings = embeddings[1]
      print(paste0("===== Computing UMAP using ",embeddings,". ====="))
      umap.genes <- AIT.anndata$var[[embeddings]]
      if(!exists("t.normalized.expr")){
        t.normalized.expr <- Matrix::t(normalized.expr)         
        t.normalized.expr <- as(t.normalized.expr,"dgCMatrix")
      }
      pcs        <- prcomp(t.normalized.expr[umap.genes,], scale = TRUE)$rotation
      embeddings <- umap(pcs[,1:number.of.pcs])$layout
      embeddings <- as.data.frame(embeddings)
      rownames(embeddings) <- rownames(counts)
    }
  }
  
  ## Add embeddings (dealing with the "X_" header)
  if(!is.null(embeddings)){
    print("===== Adding provided embeddings. =====")
    if(is.matrix(embeddings))
      embeddings <- as.data.frame(embeddings)
    if(is.list(embeddings)){
      if(is.data.frame(embeddings)){
        AIT.anndata$obsm[["X_default_standard"]] = embeddings
      } else {
      # Check the names and add "X_" if needed
      nm <- names(embeddings)
      for (i in 1:length(nm)) if(substr(nm[i],1,2)!="X_") nm[i] = paste0("X_",nm[i])
      names(embeddings) <- nm
      for(embedding in names(embeddings))
        AIT.anndata$obsm[[embedding]] = embeddings[[embedding]]
      }
    }

    ## Set default_embedding
    if(is.null(default_embedding)){
      AIT.anndata$uns$default_embedding = names(AIT.anndata$obsm)[1]
    } else {
      # Check the name and add "X_" if needed
      default_embedding = default_embedding[1]
      if(substr(default_embedding,1,2)!="X_") default_embedding = paste0("X_",default_embedding)
      default_embedding = intersect(default_embedding,names(AIT.anndata$obsm))
      if(length(default_embedding)==1){
        AIT.anndata$uns$default_embedding = default_embedding
      } else {
        AIT.anndata$uns$default_embedding = names(AIT.anndata$obsm)[1]
      }
    }
    print(paste("Default embedding set as:", substr(AIT.anndata$uns$default_embedding,1,1000)))
  }
  
  ## Add additional uns variables
  if (length(uns.variables)>0){
    for (nm in names(uns.variables))
     AIT.anndata$uns[[nm]] <- uns.variables[[nm]]
  }

  ## Write the Allen Institute Taxonomy object  # THIS SEEMS UNNECESSARY!
  print("===== Writing taxonomy anndata =====")
  AIT.anndata$write_h5ad(file.path(AIT.anndata$uns$taxonomyDir, paste0(AIT.anndata$uns$title, ".h5ad")))

  ## Read in the AIT object to confirm it was written correctly and to set parameters for the next steps  # THIS SEEMS UNNECESSARY!
  AIT.anndata = read_h5ad(file.path(AIT.anndata$uns$taxonomyDir, paste0(AIT.anndata$uns$title, ".h5ad")))
  
  ## Add dendrogram markers and membership tables, if requested
  if(add.dendrogram.markers){
    print("===== Adding dendrogram markers and membership tables for tree mapping =====")
    tryCatch({
      AIT.anndata = addDendrogramMarkers(AIT.anndata, 
                                          mode="standard", 
                                          celltypeColumn = celltypeColumn,
                                          taxonomyDir=AIT.anndata$uns$taxonomyDir,
                                          ...)
    }, error = function(e) {
      print("===== Error adding dendrogram markers. Skipping this step. =====")
      print(e)
    })
  }
  
  ## NOTE: REMOVE CHECKS OF HEIRARCHY AND DO ELSEWEHRE
  ## Add MapMyCells (hierarchical mapping) functionality, if requested
  if(addMapMyCells) {
    print("===== Adding MapMyCells (hierarchical mapping) functionality =====")
    tryCatch({
      AIT.anndata = mappingMode(AIT.anndata, mode="standard")
      AIT.anndata = addMapMyCells(AIT.anndata, hierarchy, force=TRUE)
    }, error = function(e) {
      print("===== Error adding MapMyCells functionality. Skipping this step. =====")
      print(e)
    })
  }

  ## Write the Allen Institute Taxonomy object without the normalized data (it can be recalculated on load)
  if(!is.null(AIT.anndata$X)){
    print("===== Writing taxonomy anndata without saved normalized data=====")
    X <- AIT.anndata$X
    AIT.anndata$X = NULL
    AIT.anndata$write_h5ad(file.path(AIT.anndata$uns$taxonomyDir, paste0(AIT.anndata$uns$title, ".h5ad")))
    AIT.anndata$X <- X
  }
  
  ## Check whether the taxonomy is a valid scrattch.taxonomy format
  if(check.taxonomy){
    print("===== Checking taxonomy for adherence to schema =====")
    AIT.anndata = checkTaxonomy(AIT.anndata, print.messages=print.messages)
  }
  
  ## Return the anndata object
  return(AIT.anndata)
}
