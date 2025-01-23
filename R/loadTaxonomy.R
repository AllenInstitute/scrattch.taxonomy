#' Read in a reference data set in Allen taxonomy format
#'
#' @param taxonomyDir Directory containing the Shiny taxonomy -OR- a direct h5ad file name.
#' @param anndata_file File name of the anndata object (.h5ad) to be loaded.
#' @param log.file.path Path to write log file to. Defaults to current working directory. 
#'
#' @return Organized reference object ready for mapping against.
#'
#' @export
loadTaxonomy = function(taxonomyDir, 
                        anndata_file = "AI_taxonomy.h5ad",
                        log.file.path=getwd(),
                        force=FALSE){

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

  ## Set scrattch.mapping to default standard mapping mode
  AIT.anndata$uns$mode = "standard"
  AIT.anndata$uns$taxonomyDir = taxonomyDir
  AIT.anndata$uns$title = anndata_file

  ## Return
  return(AIT.anndata)
}
