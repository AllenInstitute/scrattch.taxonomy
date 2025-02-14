#' This function builds the minimum files required for Shiny
#' Precomputed clusters must be provided.  In the anndata object these will be stored using the term "cluster".  If hierarchy[[-1]] is anything other than cluster, then any existing "cluster" column will be overwritten by hierarchy[[-1]].  Values can be provided without colors and ids (e.g., "cluster") or with them (e.g., "cluster_label" + "cluster_color" + "cluster_id").  In this case cluster_colors is ignored and colors are taken directly from the metadata.  Cluster_id's will be overwritten to match dendrogram order.
#'
#' @param meta.data Meta.data corresponding to count matrix. Rownames must be equal to colnames of counts. "clusters" must be provided (see hierarchy[[-1]] and notes).
#' @param title The file name to assign for the Taxonomy h5ad.
#' @param counts A count matrix in sparse format: dgCMatrix.
#' @param highly_variable_genes Set of features defined as highly variable genes. Provide either as a named list of vectors, or as a single vector (in which case the name "highly_variable_genes_standard" will be used).
#' @param marker_genes Set of features defined as marker genes. Provide either as a named list of vectors, or as a single vector (in which case the name "marker_genes_standard" will be used).
#' @param ensembl_id A vector of ensemble ids corresponding to the gene symbols in counts.
#' @param cluster_stats A matrix of median gene expression by cluster. Cluster names must exactly match meta.data$cluster.  If provided, will get saved to "varm$cluster_id_median_expr_[mode]"
#' @param embeddings Dimensionality reduction coordinate data.frame with 2 columns. Rownames must be equal to colnames of counts.  Either provide as a named list or as a single data.frame (in which case the name "default_standard" will be used).
#' @param dend Existing dendrogram associated with this taxonomy (e.g., one calculated elsewhere).  If NULL (default) a new dendrogram will be calculated based on the input `feature.set`
#' @param taxonomyDir The location to save Shiny objects, e.g. "/allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/10x_seq/NHP_BG_20220104/"
#' @param hierarchy List of term_set_labels in the Taxonomy ordered from most gross to most fine (e.g., subclass, supertype, neighborhood, class).
#' @param subsample The number of cells to retain per cluster (default = 2000)
#' @param cluster_colors An optional named character vector where the values correspond to colors and the names correspond to celltypes in hierarchy[[-1]].  If this vector is incomplete, a warning is thrown and it is ignored. cluster_colors can also be provided in the metadata (see notes)
#' @param default_embedding A string indicating which embedding to use for calculations.  Default (NULL) is to take the first one provided in embeddings.
#' @param uns.variables If provided, a list of additional variables to be included in the uns.  See Notes for schema variables not otherwise accounted for.
#' @param reorder.dendrogram Should dendogram attempt to match a preset order? (Default = FALSE).  If TRUE, the dendrogram attempts to match the celltype factor order as closely as possible (if celltype is a character vector rather than a factor, this will sort clusters alphabetically, which is not ideal).
#' @param add.dendrogram.markers If TRUE (default), will also add dendrogram markers to prep the taxonomy for tree mapping
#' @param addMapMyCells If TRUE (default), will also prep this taxonomy for hierarchical mapping
#' @param check.taxonomy description
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
buildTaxonomy = function(meta.data,
                         title,
                         hierarchy,
                         ## X
                         counts = NULL, ## matrix, sparse.
                         normalized.expr = NULL, ## matrix, sparse.
                         ## var
                         highly_variable_genes = NULL, ## named list
                         marker_genes = NULL, ## named list
                         ensembl_id = NULL,
                         ## varm
                         cluster_stats = NULL,
                         ## obsm
                         embeddings = NULL,
                         ## uns
                         dend = NULL,
                         taxonomyDir = getwd(),
                         cluster_colors = NULL, ## @Jeremy please handle this.
                         default_embedding = NULL,
                         uns.variables = list(),  ## @Nelson: this is my inelegant way of handling Additional uns.variables (see notes above... feel free to edit)
                         ## Additional parameters
                         subsample=2000,
                         features.dendrogram = NULL,
                         reorder.dendrogram = FALSE,
                         add.dendrogram.markers = TRUE,
                         addMapMyCells = TRUE,
                         check.taxonomy = TRUE,
                         ...){

  ## Sanity check and cleaning of parameters
  clean.params = .checkBuildTaxonomyParams(counts, 
                                            normalized.expr,
                                            meta.data, 
                                            highly_variable_genes,
                                            marker_genes,
                                            embeddings, 
                                            hierarchy,
                                            cluster_stats,
                                            taxonomyDir, 
                                            title, 
                                            dend)

  ## Pull the finest level cell type column
  celltypeColumn = names(hierarchy)[[-1]]

  ## ----------
  ## Subsample nuclei per cluster, max of subsample cells per cluster
  kpSub = subsample_taxonomy(meta.data[[celltypeColumn]], rownames(meta.data), dend, subsample)

  ## Define clusterInfo, which is used to convert cell types to subclass / neighborhood / class
  clusterInfo = meta.data %>% distinct(.data[[celltypeColumn]], .keep_all = TRUE)
  rownames(clusterInfo) = clusterInfo[[celltypeColumn]]

  ## Computing TPM matrix if count matrix exists
  if(!is.null(counts) & is.null(normalized.expr)){
    print("===== Normalizing count matrix to log2(CPM) =====")
    normalized.expr = as(scrattch.bigcat::logCPM(counts), "dgCMatrix")
  }else{
    print("===== No provided count matrix. Skipping TPM calculation. =====")
  }

  ## Compute cluster medians if needed
  if(!is.null(normalized.expr)){
    if(is.null(cluster_stats)){
      print("===== Computing median expr. at taxonomy leaves =====")
      ## Get cluster medians
      cluster   = meta.data[[celltypeColumn]]; names(cluster) = rownames(meta.data)
      cluster_stats = scrattch.bigcat::get_cl_medians(normalized.expr, cluster)
    }
  }

  ##  Build the dendrogram
  print("===== Building dendrogram =====")
  if(!is.null(schema$dendrogram)){
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
    invisible(capture.output({  # Avoid printing lots of numbers to the screen
      dend.result = build_dend(
        cl.dat  = cluster_stats[feature.set,],
        cl.cor  = NULL,
        # l.color = use.color, @Jeremy please handle this.
        l.rank  = l.rank, 
        nboot   = 1,
        ncores  = 1)
    }))
    print("...dendrogram built.")
  }

  ## If highly_variable_genes and/or marker genes are provided as vectors (as with previous versions of scrattch.taxonomy, convert to named lists)
  if(is.character(highly_variable_genes)){
    highly_variable_genes <- list(highly_variable_genes_standard = highly_variable_genes)
  }
  if(is.character(marker_genes)){
    marker_genes <- list(marker_genes_standard = marker_genes)
  }
  
  ## Compute most likely "ontology_term_id" for "organism", "anatomical_region", "self_reported_sex", "self_reported_ethnicity", "assay", and "disease"
  # meta.data <- computeOntologyTerms(meta.data)
  # NOTE this is a wrapper function for updateTaxonomyMetadata and found in updateTaxonomyMetadata.R
  # I think it's better that this gets called outside this function
  
  ## Build the AIT object
  print("===== Building taxonomy anndata =====")
  AIT.anndata = AnnData(
    X = if(!is.null(normalized.expr)) Matrix::t(normalized.expr) else NULL, ## logCPM ensured to be in sparse column format
    raw = if(!is.null(counts)) list(X = as(Matrix::t(counts), "dgCMatrix"), var = data.frame("gene" = rownames(counts))) else NULL, ## Store counts matrix
    obs = meta.data,
    var = if(!is.null(counts))
                    data.frame("gene" = rownames(counts), 
                               row.names = rownames(counts))
          else data.frame(),
    varm = list(
      "cluster_id_median_expr_standard" = cluster_stats ## Median expression by cluster
    ),
    obsm = list(  ## A data frame with cell_id, and 2D coordinates for umap (or comparable) representation(s)
    ),  
    uns = list(
      dend       = list("standard" = toJSON(dend_to_json(dend.result$dend))), ## JSON dendrogram
      filter     = list("standard" = !(colnames(counts) %in% kpSub)), ## Filtered cells
      mapmycells = list("standard" = list()),
      mode = "standard", ## Default mode to standard
      cellSet = rownames(meta.data),
      clusterInfo = clusterInfo,
      clusterStatsColumns = list("standard" = colnames(cluster_stats)),
      title = title,
      hierarchy = hierarchy,
      taxonomyDir = file.path(normalizePath(taxonomyDir), leading_string="/"), ## Normalize path in case where user doesn't provide absolute path.
      schema_version = as.character(packageVersion("scrattch.taxonomy"))  # Tying schema version to scrattch.taxonomy version. @Nelson, I'm open to another option.
    )
  )

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

  ## Add ensemble_id into AIT object
  if(!is.null(ensembl_id)){
    if(dim(AIT.anndata$var)[1]==length(ensembl_id)){
      AIT.anndata$var$ensembl_id = as.character(ensembl_id)
    } else {
      print("===== Provided ensembl_id vector of different length from gene length in var. Skipping adding ensembl_ids. =====")
    }
  }
  
  ## Add embeddings (dealing with the "X_" header)
  if(!is.null(embeddings)){
    print("===== Adding provided embeddings. =====")
    if(is.list(embeddings)){
      # Check the names and add "X_" if needed
      nm <- names(embeddings)
      for (i in 1:length(nm)) if(substr(nm[i],1,2)!="X_") nm[i] = paste0("X_",nm[i])
      names(embeddings) <- nm
      for(embedding in names(embeddings))
        AIT.anndata$obsm[[embedding]] = embeddings[[embedding]]
    } else {
      if(is.data.frame(embeddings))
        AIT.anndata$obsm[["X_default_standard"]] = embeddings
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
    print(paste("Default embedding set as:",substr(AIT.anndata$uns$default_embedding,3,1000)))
  }
  
  ## Add additional uns variables
  if (length(uns.variables)>0){
    for (nm in names(uns.variables))
     AIT.anndata$uns[[nm]] <- uns.variables[[nm]]
  }

  ## Write the Allen Institute Taxonomy object
  print("===== Writing taxonomy anndata =====")
  AIT.anndata$write_h5ad(file.path(AIT.anndata$uns$taxonomyDir, paste0(AIT.anndata$uns$title, ".h5ad")))

  ## Read in the AIT object to confirm it was written correctly and to set parameters for the next steps
  AIT.anndata = read_h5ad(file.path(AIT.anndata$uns$taxonomyDir, paste0(AIT.anndata$uns$title, ".h5ad")))
  
  ## Add dendrogram markers and membership tables, if requested
  if(add.dendrogram.markers){
    print("===== Adding dendrogram markers and membership tables for tree mapping =====")
    if(is.character(hierarchy)) hierarchy <- as.list(hierarchy)
    AIT.anndata = addDendrogramMarkers(AIT.anndata, 
                                        mode="standard", 
                                        celltypeColumn = celltypeColumn,
                                        taxonomyDir=AIT.anndata$uns$taxonomyDir,
                                        ...)
    # The reference probability matrix for the subsetted taxonomy is defined and outputted in this function as well
    # $memb[[mode.name]]
    # ...$memb.ref,
    # ...$map.df.ref
  }
  
  ## NOTE: REMOVE CHECKS OF HEIRARCHY AND DO ELSEWEHRE
  ## Add MapMyCells (hierarchical mapping) functionality, if requested
  if(addMapMyCells) {
    if(is.character(hierarchy)) hierarchy <- as.list(hierarchy)
    if((length(hierarchy)==0)|(sum(class(hierarchy)=="list")<1)){
      warning("hierarchy must be a list of term_set_labels in the reference taxonomy ordered from most gross to most fine included in AIT_anndata or provided separately. Since this is NOT the case, addMapMyCells is being skipped")
    } else{
      print("===== Adding MapMyCells (hierarchical mapping) functionality =====")
      AIT.anndata = mappingMode(AIT.anndata, mode="standard")
      AIT.anndata = addMapMyCells(AIT.anndata, names(hierarchy), force=TRUE)
    }
  }

  ## Check whether the taxonomy is a valid scrattch.taxonomy format
  AIT.anndata$uns$check = list(isValid=NA)
  if(check.taxonomy){
    AIT.anndata$uns$check <- checkTaxonomy(AIT.anndata)
  }
  
  ## Return the anndata object
  return(AIT.anndata)
}
