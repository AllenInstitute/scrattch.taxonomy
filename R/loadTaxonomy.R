#' Read in a reference data set in Allen taxonomy format
#'
#' @param taxonomyDir Directory containing the AIT file -OR- a direct h5ad file name -OR- a URL of a publicly accessible AIT file.
#' @param anndata_file If taxonomyDir is a directory, anndata_file must be the file name of the anndata object (.h5ad) to be loaded in that directory. If taxonomyDir is a file name or a URL, then anndata_file is ignored.
#' @param log.file.path Path to write log file to. Defaults to current working directory. 
#'
#' @return Organized reference object ready for mapping against.
#' 
#' @import anndata
#'
#' @export
loadTaxonomy = function(taxonomyDir = getwd(), 
                        anndata_file = "AI_taxonomy.h5ad",
                        log.file.path=getwd(),
                        force=FALSE){
  
  ## Allow for h5ad as the first/only input
  if(grepl("h5ad", taxonomyDir)){
    anndata_file = taxonomyDir
    taxonomyDir  = getwd()
  }
  ## Make sure the taxonomy path is an absolute path
  taxonomyDir = normalizePath(taxonomyDir, winslash = "/")
  
  ## If anndata_file is a URL, (1) parse the bucket out, (2) check whether the file is currently in the working directory, 
  ##   (3) download it if not, and then (4) rename anndata_file to the file name.
  # https://allen-cell-type-taxonomies.s3.us-west-2.amazonaws.com/Mouse_VISp_ALM_SMART_seq_04042025.h5ad
  if(grepl("http", anndata_file)&grepl("s3",anndata_file)){
    anndata_parts  <- strsplit(anndata_file, "/")[[1]]
    anndata_object <- anndata_parts[length(anndata_parts)]
    if(!file.exists(anndata_object)){
      download.file(anndata_file, file.path(taxonomyDir, anndata_object), timeout=9000)
    }
    anndata_file = anndata_object
  }
  
  ## Load from directory name input 
  if(file.exists(file.path(taxonomyDir, anndata_file))){
    print("Loading reference taxonomy into memory from .h5ad")
    ## Load taxonomy directly!
    AIT.anndata = read_h5ad(file.path(taxonomyDir, anndata_file))
    ## Default mode is always standard
    AIT.anndata$uns$mode = "standard"
    ##
    # for(mode in names(AIT.anndata$uns$dend)){
    #   invisible(capture.output({
    #     if(grepl("dend.RData", AIT.anndata$uns$dend[[mode]])){
    #       print("Loading an older AIT .h5ad version. Converting dendrogram to JSON format for mapping.")
    #       dend = readRDS(AIT.anndata$uns$dend[[mode]])
    #       AIT.anndata$uns$dend[[mode]] = toJSON(dend_to_json(dend))
    #     }
    #   }))
    # }
    if(("taxonomyName" %in% colnames(AIT.anndata$obs)) & (!"title" %in% colnames(AIT.anndata$obs))){
      AIT.anndata$obs$title = anndata_file$obs$taxonomyName
    }
    ## Ensure anndata is in scrattch.mapping format
    # if(!checkTaxonomy(AIT.anndata,log.file.path)){
    #  stop(paste("Taxonomy has some breaking issues.  Please check checkTaxonomy_log.txt in", log.file.path, "for details"))
    # }
  }else{
    stop("Required files to load Allen Institute taxonomy are missing.")
  }
  
  ## If counts are included but normalized counts are not, calculate normalized counts
  if((!is.null(AIT.anndata$raw$X))&(is.null(AIT.anndata$X))){
    normalized.expr = log2CPM_byRow(AIT.anndata$raw$X)
    AIT.anndata$X   = normalized.expr 
  }

  ## Set scrattch.mapping to default standard mapping mode
  AIT.anndata$uns$mode = "standard"
  AIT.anndata$uns$taxonomyDir = taxonomyDir
  AIT.anndata$uns$title = anndata_file

  ## Return
  return(AIT.anndata)
}
