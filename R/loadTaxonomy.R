#' Read in a reference data set in Allen taxonomy format
#'
#' @param taxonomyDir Directory containing the Shiny taxonomy -OR- a direct h5ad file name.
#' @param anndata_file File name of the anndata object (.h5ad) to be loaded.
#' @param sample_id Field in reference taxonomy that defines the sample_id.
#' @param hGenes User supplied variable gene vector.  If not provided, then all genes are used.
#' @param gene_id Field in counts.feather that defines the gene_id.
#' @param log.file.path Path to write log file to. Defaults to current working directory. 
#' @param force Force rebuild the anndata object for the taxonomy.
#'
#' @return Organized reference object ready for mapping against.
#'
#' @export
loadTaxonomy = function(taxonomyDir, 
                        anndata_file = "AI_taxonomy.h5ad",
                        sample_id = "sample_id", 
                        hGenes=NULL, 
                        gene_id = "gene",
                        log.file.path=getwd(),
                        force=FALSE){

  ## Load from directory name input 
  if(file.exists(file.path(taxonomyDir, anndata_file)) & force == FALSE){
    print("Loading reference taxonomy into memory from .h5ad")
    ## Load taxonomy directly!
    AIT.anndata = read_h5ad(file.path(taxonomyDir, anndata_file))
    ## Default mode is always standard
    AIT.anndata$uns$mode = "standard"
    ##
    for(mode in names(AIT.anndata$uns$dend)){
      invisible(capture.output({
        if(grepl("dend.RData", AIT.anndata$uns$dend[[mode]])){
          print("Loading an older AIT .h5ad version. Converting dendrogram to JSON format for mapping.")
          dend = readRDS(AIT.anndata$uns$dend[[mode]])
          AIT.anndata$uns$dend[[mode]] = toJSON(dend_to_json(dend))
        }
      }))
    }
    ## Ensure anndata is in scrattch.mapping format
    if(!checkTaxonomy(AIT.anndata,log.file.path)){
     stop(paste("Taxonomy has some breaking issues.  Please check checkTaxonomy_log.txt in", log.file.path, "for details"))
    }
  } else if(all(file.exists(c(file.path(taxonomyDir,"anno.feather"), 
                              file.path(taxonomyDir,"data.feather"), 
                              file.path(taxonomyDir,"counts.feather"), 
                              file.path(taxonomyDir,"tsne.feather"))))){
    ##
    print("Loading reference taxonomy into memory from .feather")
    
    ## Read in reference data and annotation files and format correctly
    annoReference   = feather(file.path(taxonomyDir,"anno.feather")) 
    exprReference   = feather(file.path(taxonomyDir,"data.feather"))
    
    ## Convert log2CPM-normalized data into a matrix
    datReference = as.matrix(exprReference[,names(exprReference)!=sample_id])

    ## Read in reference count data if available
    if(file.exists(file.path(taxonomyDir,"counts.feather"))){
      print("Loading counts matrix")
      countsReference = feather(file.path(taxonomyDir,"counts.feather"))  # Note that this is transposed relative to data.feather
      ## Convert count data into a matrix
      countsReference = as.matrix(countsReference[,names(countsReference)!=gene_id])
      rownames(countsReference) = colnames(datReference)
      countsReference = Matrix::t(countsReference) ## Put into cell x gene format
    }else{
      countsReference = NULL
    }
    
    ## Match meta.data to data
    annoReference = as.data.frame(annoReference[match(exprReference[[sample_id]], annoReference[[sample_id]]),])
    rownames(annoReference) = rownames(datReference) = annoReference[[sample_id]]
    
    ## Consider only genes present in both data sets
    if(!is.null(hGenes)){ 
      feature.set = intersect(hGenes, colnames(datReference)) 
    } else { 
      feature.set = colnames(datReference) 
    }
    
    ## Read in cluster info
    clustersUse = unique(annoReference$cluster_label)
    clusterInfo = as.data.frame(annoReference) ## No dendrogram ordering so just convert to df
 
    ## Read in the umap
    if(!file.exists(file.path(taxonomyDir,"tsne.feather"))){
      umap.coords = NULL
    }else{
      tryCatch(
        expr = {
          umap.coords = as.data.frame(read_feather(file.path(taxonomyDir,"tsne.feather")))
          rownames(umap.coords) = umap.coords[,sample_id]
          umap.coords = umap.coords[rownames(annoReference),]
        },
        error = function(e){ 
          print("Error caught for umap loading. Setting all UMAP values to 0.")
          l = length(rownames(annoReference))
          umap.coords = data.frame(sample_id=rownames(annoReference),all_x=rep(0,l),all_y=rep(0,l))
          rownames(umap.coords) = umap.coords[,sample_id]
        },
        warning = function(w){
        },
        finally = {
          print("Done loading UMAP")
        }
      )
    }
    
    # This next bit should NOT be needed, but the anndata crashes without it
    if(is_tibble(umap.coords)){
      umap.coords = as.data.frame(umap.coords)
      rownames(umap.coords) = umap.coords[,sample_id]
    }

    ##
    dend = readRDS(file.path(taxonomyDir, "dend.RData", leading_string="/"))

    ##
    AIT.anndata = AnnData(
      X = datReference, ## logCPM
      obs = annoReference,
      var = data.frame("gene" = colnames(datReference), 
                       "highly_variable_genes" = colnames(datReference) %in% feature.set, 
                       row.names=colnames(datReference)),
      layers = list(
        counts = countsReference ## Count matrix
      ),
      obsm = list(
        umap = umap.coords ## A data frame with sample_id, and 2D coordinates for umap (or comparable) representation(s)
      ),
      uns = list(
        dend        = list("standard" = toJSON(dend_to_json(dend))), # FILE NAME with dendrogram
        filter      = list("standard" = rep(FALSE, nrow(datReference))),
        QC_markers  = list("standard" = NA), ## Standard should always be empty for QC_markers
        mode = "standard", ## Default mode to standard
        clustersUse = clustersUse,
        clusterInfo = clusterInfo,
        taxonomyName = gsub(".h5ad","", anndata_file),
        taxonomyDir = file.path(taxonomyDir, leading_string="/")
      )
    )
  }else{
    stop("Required files to load Allen Institute taxonomy are missing.")
  }

  ## Set scrattch.mapping to default standard mapping mode
  AIT.anndata$uns$mode = "standard"

  ## Return
  return(AIT.anndata)
}
