#' Create a taxonomy mode for mapping to taxonomy subset (or child)
#'
#' This function creates a new mode for mapping to a subset of a taxonomy. This includes filtering cells, setting new variable genes and/or marker genes, calculating mapping statistics, 
#'
#' @param AIT.anndata A reference taxonomy anndata object.
#' @param mode.name A name to identify the new taxonomy version.
#' @param retain.cells A boolean vector of length [number of cells] indicating which cells should be retained (TRUE) or filtered (FALSE) -OR- a character vector with sample names indicating which cells should be retained. Default is to retain the cells included in "stadard" mode.
#' @param retain.clusters A character vector with cluster names (e.g., values in the "cluster_id" column) indicating which clusters should be retained. Default is to retain all clusters with at least 2 retained cells (clusters with exactly 1 cell can cause some functions to crash).
#' @param subsample The number of cells to retain per cluster (default = 100). Note that subsampling happens AFTER retain.cells and retail.clusters filtering.
#' @param highly_variable_genes Set of features defined as highly variable genes OR a number of binary genes to calculate (we recommend ~1000 - ~5000, for <100 to ~5000 cell types). If a feature list is provided, provide either as a named list of vectors, or as a single vector (in which case the name "highly_variable_genes_[mode.name]" will be used). "highly_variable_genes_[mode.name]" will also be used for calculated variable genes. Optional input, but for proper mapping we recommend including either highly_variable_genes or marker_genes. If nothing is provided (default=NULL), standard mode markers will be used for mapping algorithms based on these gene lists and may cause problems.
#' @param marker_genes Set of features defined as marker genes. Provide either as a named list of vectors, or as a single vector (in which case the name "marker_genes_[mode.name]" will be used).
#' @param embeddings Dimensionality reduction coordinate data.frame with 2 columns or a string with the column name for marker_genes or variable_genes from which a UMAP should be calculated. If coordinates are provided, rownames must be equal to colnames of counts.  Either provide as a named list or as a single data.frame (in which case the name "default_[mode.name]" will be used). Optional - if nothing is provided (default=NULL) the relevant subset of the default standard embedding will be used.
#' @param number.of.pcs Number of principle components to use for calculating UMAP coordinates (default=30). This is only used in embeddings corresponds to a variable gene column from which a UMAP should be calculated.
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
                             number.of.pcs = 30,
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
  AIT.anndata$uns$clusterStatsColumns[[mode.name]] = mode.clusters

  
  ## Filter stats files if they exist
  #AIT.anndata$varm[[paste0("cluster_id_median_expr_",mode.name)]] =  AIT.anndata$varm[["cluster_id_median_expr_standard"]][,mode.clusters]
  #  THIS IS REDUNDANT TO SAVE, SO WE WON'T.  IT CAN BE SUBSET AS NEEDED
  
  
  ## DEAL WITH MODE-SPECIFIC HIGHLY VARIABLE GENES
  
  ## highly_variable_genes is a data.frame with gene names in rows and various sets in columns
  if(!is.null(highly_variable_genes)){
    # If numeric, calculating the top N binary genes and save as a list first
    if(is.numeric(highly_variable_genes)){
      counts <- AIT.anndata$raw$X
      highly_variable_genes <- round(max(min(dim(counts)[1],highly_variable_genes[1]),100))  # Make sure the numeric value is legal
      counts <- BiocGenerics::t(counts[retain,]) # Transpose again.  NOT IDEAL!
      counts <- as(counts, "dgCMatrix")
      rownames(counts) <- rownames(AIT.anndata$var)
      colnames(counts) <- rownames(AIT.anndata$obs)[retain]
      binary.genes          <- top_binary_genes(counts, as.character(AIT.anndata$obs[[celltypeColumn]][retain]), highly_variable_genes)
      highly_variable_genes <- list(binary.genes)
      names(highly_variable_genes) <- paste0("highly_variable_genes_",mode.name)
    }
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
  
  
  ## MODIFY THE DENDROGRAM AND SAVE IN THE ANNDATA
  ### -- FUTURE UPDATE: Allow a separate dendrogram to be entered as a variable rather than automatically subsetting the existing dendrogram
  ### -- FUTURE UPDATE: Allow for non-binary dendrograms
  
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
  
  
  ## DEAL WITH MODE-SPECIFIC EMBEDDINGS
  
  ## Calculate embeddings if requested
  if(is.character(embeddings)){
    embeddings <- intersect(embeddings,colnames(AIT.anndata$var))
    if(is.null(embeddings)){
      normalized.expr <- NULL
    } else {
      normalized.expr <- as(scrattch.bigcat::logCPM(counts), "dgCMatrix")
    }
    if((!is.null(embeddings))&(!is.null(normalized.expr))){
      embeddings = embeddings[1]
      print(paste0("===== Computing UMAP using ",embeddings,". ====="))
      umap.genes <- AIT.anndata$var[[embeddings]]
      pcs        <- prcomp(normalized.expr[umap.genes,], scale = TRUE)$rotation
      embeddings <- umap(pcs[,1:number.of.pcs])$layout
      embeddings <- as.data.frame(embeddings)
      rownames(embeddings) <- colnames(counts)
      
      # Add 0s in non-retained cells so it will fit in the obsm slot
      embeddings2 <- matrix(0,nrow=length(rownames(AIT.anndata$obs)),ncol=2)
      rownames(embeddings2) <- rownames(AIT.anndata$obs)
      colnames(embeddings2) <- colnames(embeddings2)
      embeddings2[colnames(counts),]  <- as.matrix(embeddings)
      embeddings <- embeddings2
    }
  }
  
  ## Add embeddings (dealing with the "X_" header)
  if(!is.null(embeddings)){
    print("===== Adding provided embeddings. =====")
    if(is.matrix(embeddings))
      embeddings <- as.data.frame(embeddings)
    if(is.list(embeddings)){
      if(is.data.frame(embeddings)){
        AIT.anndata$obsm[[paste0("X_default_",mode.name)]] = embeddings
      } else {
        # Check the names and add "X_" if needed
        nm <- names(embeddings)
        for (i in 1:length(nm)) if(substr(nm[i],1,2)!="X_") nm[i] = paste0("X_",nm[i])
        names(embeddings) <- nm
        for(embedding in names(embeddings))
          AIT.anndata$obsm[[embedding]] = embeddings[[embedding]]
      }
    }
    
  }
   
   
  ## MAPMYCELLS FUNCTIONALITY
  
  ## Add MapMyCells (hierarchical mapping) functionality, if requested
  if(addMapMyCells) {
    print("===== Adding MapMyCells (hierarchical mapping) functionality =====")
    tryCatch({
      current.mode = AIT.anndata$uns$mode
      AIT.anndata  = mappingMode(AIT.anndata, mode=mode.name)
      AIT.anndata  = addMapMyCells(AIT.anndata, hierarchy, force=TRUE)
      AIT.anndata  = mappingMode(AIT.anndata, mode=current.mode)
    }, error = function(e) {
      print("===== Error adding MapMyCells functionality. Skipping this step. =====")
      print(e)
    })
  }
  
  ## Write the Allen Institute Taxonomy object without the normalized data (it can be recalculated on load)
  print("===== Writing taxonomy anndata without saved normalized data=====")
  AIT.anndata2 = AIT.anndata
  AIT.anndata2$X = NULL
  AIT.anndata2$uns$title <- gsub(".h5ad","",AIT.anndata2$uns$title)
  AIT.anndata2$write_h5ad(file.path(AIT.anndata2$uns$taxonomyDir, paste0(AIT.anndata2$uns$title, ".h5ad")))
  rm(AIT.anndata2)
  
  ##
  return(AIT.anndata)

}
