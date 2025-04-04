#' Function to define AIT schema
#'
#' @return AIT schema
#'
#' @keywords internal
.schemaAIT = function(){
  data = data.frame("")
}


#' Convert to sentence case
#'
#' @param str A character string
#'
#' @return The character string in sentence case
tosentence <- function(str) {
  str <- as.character(str)
  paste0(toupper(substr(str,1,1)),tolower(substr(str,2,nchar(str))))
}



#' Convert the first character of a string to lowercase
#'
#' @param str A character string
#'
#' @return The character string with specific characters changed to lower case
firsttolower <- function(str) {
  str <- as.character(str)
  paste0(tolower(substr(str,1,1)),substr(str,2,nchar(str)))
}


#' Convert gene symbols to Ensembl IDs
#' 
#' NOTE: This function requires an internet connection.
#' 
#' @param gene.symbols List of gene symbols to convert to Ensembl IDs.
#' @param ncbi.taxid The integer part of the NCBITaxon ID for the species want to convert genes between.
#' @param use.synonyms If TRUE (default) will search synonyms for current gene symbols to try and match Ensembl IDs
#' @param remove.duplicates If TRUE (default) any genes that share Ensembl IDs with any other genes will have their Ensembl IDs set to NA to avoid ambiguity. In cases where a duplicate is introduced in the synonyms, the original gene.symbol will retain the Ensembl IDs and synonym duplicates will be set to NA.
#' @param includeNonMammalianSpecies Default (FALSE) only considers mammalian species. Set to TRUE if non-mammalian species are considered (much slower).
#'
#' @import data.table
#'
#' @return Ensembl IDs for the input gene.symbols (or NA if not found or if duplicated)
#'
#' @export
geneSymbolToEnsembl <- function (gene.symbols, ncbi.taxid = 9606, use.synonyms = TRUE, remove.duplicates=TRUE, includeNonMammalianSpecies=FALSE)
{
  ## Read the gene information file from the smallest available file
  if(ncbi.taxid==9606){  # If we go back to reporting for only primary species, we can speed up the code by putting this back
    geneInfo <- fread("https://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/Mammalia/Homo_sapiens.gene_info.gz")
  } else if(ncbi.taxid==10090){
    geneInfo <- fread("https://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/Mammalia/Mus_musculus.gene_info.gz")
  } else if (!includeNonMammalianSpecies){
    geneInfo <- fread("https://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/Mammalia/All_Mammalia.gene_info.gz")
    geneInfo <- geneInfo[geneInfo$`#tax_id`==ncbi.taxid,]
  } else {
    print("Non-mammalian species indicated, so this step may be very slow...")
    geneInfo <- fread("https://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/All_Data.gene_info.gz") # a 1GB file!
    geneInfo <- geneInfo[geneInfo$`#tax_id`==ncbi.taxid,]
  } 
  if(dim(geneInfo)[1]==0) stop(paste0("ERROR: No gene information available for NCBI taxid",ncbi.taxid,". Exiting."))
  
  # Extract and return ensembl information from dbXrefs 
  ensembl.all  <- as.character(lapply(geneInfo$dbXrefs, function(x) strsplit(strsplit(x,"Ensembl:")[[1]][2],"\\|")[[1]][1]))
  ensembl.true <- setNames(ensembl.all[match(gene.symbols,geneInfo$Symbol)],gene.symbols)
  
  # Return what it finds if we don't want to use synonyms
  if(!use.synonyms) return(ensembl.true)
  
  # Explore synonyms for values that are NA, if requested
  ensembl  <- ensembl.true
  synonyms <- geneInfo$Synonyms
  
  #ens <- as.character(lapply(gene.symbols[is.na(ensembl.true)], function(x) ensembl.all[grep(x,synonyms)[1]])) # This line is very slow
  
  # Start faster replacement for above line
  syn <- apply(cbind(ensembl.all,synonyms),1,function(x) { 
    a <- as.character(x)
    if(sum(is.na(x))>0) return(NULL)
    y = strsplit(x[2],"\\|")[[1]];
    as.data.frame(cbind(rep(x[1],length(y)),y))
  })
  syn2  <- list_rbind(syn)
  merge <- setNames(syn2[,1],syn2[,2])
  ens   <- merge[gene.symbols[is.na(ensembl.true)]]
  # End faster replacement for above line
  
  ensembl[is.na(ensembl.true)] <- ens
  if(!remove.duplicates) return(ensembl)
  
  # Remove duplicates introduced by synonyms, if requested
  ensembl[is.element(ensembl,names(table(ensembl)[table(ensembl)>1]))] = NA  # Remove all duplicates
  ensembl[!is.na(ensembl.true)] <- ensembl.true[!is.na(ensembl.true)]        # Return gene symbol values
  ensembl[is.element(ensembl,names(table(ensembl)[table(ensembl)>1]))] = NA  # Remove duplicates in gene symbol values
  ensembl
}



#' Determine the NCBITaxon ID for a species from the scientific name
#' 
#' This function returns the NCBITaxon ID for a species from the scientific name.  It's probably easier just looking it up via a google search, and most of the ones used at the Allen Institute are listed here: https://github.com/AllenInstitute/GeneOrthology.
#'
#' @param species Scientific name for one or more species. Not case sensitive, but make sure you spell it correctly or the function will return NA.
#' @param ncbitaxon_obo Variable in ontologyIndex format read in from obo.file file.
#' @param obo.file Location to look for or put the ncbitaxon.obo file, if ncbitaxon_obo is not provided
#' @param obo.url URL to download obo.file from, if not already downloaded. We recommend keeping the default and only switching to "http://purl.obolibrary.org/obo/ncbitaxon.obo" if absolutely necessary.
#'
#' @import ontologyIndex
#'
#' @return an AIT.anndata object with the updated/additional vector of highly variable genes
#'
#' @export
getNCBITaxon <- function(species,
                         ncbitaxon_obo=NULL,
                         obo.file="ncbitaxon.obo",
                         obo.url="https://raw.githubusercontent.com/obophenotype/ncbitaxon/refs/heads/master/subsets/taxslim.obo"
                                 # Larger option is "http://purl.obolibrary.org/obo/ncbitaxon.obo"
                         ){
  if(is.null(ncbitaxon_obo)){
    if(!file.exists(obo.file)){
      file   <- try(download.file(obo.url,obo.file))
      if("try-error" %in% class(file)) stop(paste(obo.url,"currently inexcessible and",obo.file,"does not exist.  Exiting."))
    }
    ncbitaxon_obo <- ontologyIndex::get_OBO(obo.file)
  }
  species <- tosentence(species)
  have.species <- intersect(ncbitaxon_obo$name,species)
  taxid   <- setNames(names(ncbitaxon_obo$name[match(have.species,ncbitaxon_obo$name)]),have.species)
  if(length(taxid)==0) return(setNames(rep(NA,length(species)),species))
  taxid   <- as.numeric(gsub("NCBITaxon:","",taxid))
  taxid   <- setNames(taxid,have.species)
  setNames(taxid[match(species,names(taxid))],species)  # Return in same order, with NA for missing values
}


#' Function to add or update highly variable genes for the current mode
#'
#' @param AIT.anndata A reference taxonomy anndata object.
#' @param variable.genes Set of variable genes to add to the list.name slot in obs. Can be provided as a logical or character vector.
#' @param list.name Which slot in obs to should the highly variable genes go?  Default is highly_variable_genes_[mode]
#' @param default.list.name Which slot should highly variable genes be copied from if none are provided? Default is highly_variable_genes_standard
#'
#' @return an AIT.anndata object with the updated/additional vector of highly variable genes
#'
#' @export
updateHighlyVariableGenes = function(AIT.anndata,
                                     variable.genes=NULL,
                                     mode = AIT.anndata$uns$mode,
                                     list.name = paste0("highly_variable_genes_",AIT.anndata$uns$mode),
                                     default.list.name = "highly_variable_genes_standard"){
  
  if(is.null(variable.genes)){
    if(is.null(AIT.anndata$var[,default.list.name])){
       print(paste0("updateHighlyVariableGenes is returning the starting object since AIT.anndata$var$",default.list.name," is unspecified."))
     } else {
       print("Setting the mode-specific variable genes as the global variable genes since variable.genes=NULL.")
       AIT.anndata$var[,list.name] = AIT.anndata$var[,default.list.name]
     }
    return(AIT.anndata)
  } else if (is.logical(variable.genes)) {
    if (length(variable.genes)!=dim(AIT.anndata)[2]) {
      warning("If variable.genes is logical it must be the same length as the total number of genes in AIT.anndata. updateHighlyVariableGenes is returning the starting object.")
      return(AIT.anndata)
    }
    variable.genes.vector = variable.genes
  } else if (is.character(variable.genes)){
    variable.genes <- intersect(variable.genes,AIT.anndata$var_names)
    if (length(variable.genes)<=2){
      warning("More than 2 valid gene names must be provided in variable.genes. updateHighlyVariableGenes is returning the starting object.")
      return(AIT.anndata)
    } 
    variable.genes.vector <- is.element(AIT.anndata$var_names,variable.genes)
  } else {
    warning("variable.genes must be a character or logical vector. updateHighlyVariableGenes is returning the starting object.")
    return(AIT.anndata)
  }
  AIT.anndata$var[,list.name] = variable.genes.vector
  AIT.anndata
}


#' Function to add or update marker genes for the current mode
#'
#' @param AIT.anndata A reference taxonomy anndata object.
#' @param variable.genes Set of variable genes to add to the list.name slot in obs.
#' @param list.name Which slot in obs to should the highly variable genes go?  Default is marker_genes_[mode]
#' @param default.list.name Which slot should highly variable genes be copied from if none are provided? Default is marker_genes_standard
#'
#' @return an AIT.anndata object with the updated/additional vector of marker genes
#'
#' Note that this is a wrapper of updateHighlyVariableGenes with different default names
#'
#' @export
updateMarkerGenes = function(AIT.anndata,
                             marker.genes=NULL,
                             list.name = paste0("marker_genes_",AIT.anndata$uns$mode),
                             default.list.name = "marker_genes_standard"){
  updateHighlyVariableGenes(AIT.anndata=AIT.anndata, variable.genes=marker.genes, list.name=list.name, default.list.name=default.list.name)
}
    
                         
#' Function to update meta.data
#'
#' @param meta.data A data.frame with cell metadata
#' @param cluster_colors A named vector of colors for each cluster
#' 
#' @return auto_annotated meta.data
#'
#' @keywords internal
.formatMetadata = function(meta.data, cluster_colors=NULL){

  print("===== Format metadata table for AIT schema =====")
  ## Check for duplicate columns (e.g., XXXX_label and XXXX together will crash auto_annotate)  # NEW #
  column_names = colnames(meta.data)
  column_names_revised = gsub("_label$","", column_names)
  duplicates <- names(table(column_names_revised))[table(column_names_revised)>1]
  if(length(duplicates > 0)){
    warning("Duplicate entries for", paste(duplicates, collapse=" and "),"have been DELETED. Please check output carefully!")
    meta.data <- meta.data[,setdiff(column_names, duplicates)]
  }

  ## Run auto_annotate
  meta.data = auto_annotate(meta.data, "cell_id")

  ## Convert chars and factors to characters (moved to AFTER auto_annotate, so numbers can be assigned in correct order for factors)
  for (col in colnames(meta.data)){ 
      if(is.character(meta.data[,col]) | is.factor(meta.data[,col])){
          meta.data[,col] = as.character(meta.data[,col])
      }
  }
  # Shouldn't be needed as factors are okay, and actually suggested for cluster order
  
  ## Varibow color set is broken -- this will fix it
  for (col in which(grepl("_color",colnames(meta.data)))){
      kp = nchar(meta.data[,col])==5
      meta.data[kp,col] = paste0(meta.data[kp,col],"FF")
  }

  ## Adjust the cluster colors to match cluster_colors, if available. 
  if(!is.null(cluster_colors)){
    if(length(setdiff(meta.data$cluster,names(cluster_colors)))>0){
      warning("cluster_colors is not a named vector with colors for every cluster and will therefore be ignored.")
    } else {
      meta.data$cluster_color <- as.character(cluster_colors[meta.data$cluster_label])
    }
  }
  return(meta.data)
}

#' Function to sanity check buildTaxonomy parameters
#'
#' @param counts A count matrix (cells x genes)
#' @param meta.data A data.frame with cell metadata
#' @param celltypeColumn The column name in meta.data that contains the cell type information
#' @param umap.coords A matrix of UMAP coordinates
#' 
#' @return Stops the function if any of the parameters are not as expected
#'
#' @keywords internal
.checkBuildTaxonomyParams <- function(counts, 
                                      normalized.expr,
                                      meta.data, 
                                      highly_variable_genes,  # This is checked within buildTaxonomy now
                                      marker_genes,           # This is checked within buildTaxonomy now
                                      embeddings, 
                                      celltypeColumn,
                                      cluster_stats,
                                      taxonomyDir, 
                                      title, 
                                      dend){
  if(sum(is.element(celltypeColumn, colnames(meta.data)))==0){stop("cluster column must be defined in the meta.data object")}
  if(!all(rownames(counts) == rownames(meta.data))){stop("Colnames of `counts` and rownames of `meta.data` do not match.")}
  if(!is.data.frame(meta.data)){stop("meta.data must be a data.frame, convert using as.data.frame(meta.data)")}

  ## Sanity checks on data matrices
  if(!is.null(counts)){
    if(!is(counts, 'sparseMatrix')){stop("`counts` must be a sparse matrix.")}
    if(!is.null(normalized.expr)){
      if(!is(normalized.expr, 'sparseMatrix')){stop("`normalized.expr` must be a sparse matrix.")}
      if(!all(colnames(counts) == colnames(normalized.expr))){stop("Rownames of `counts` and `normalized.expr` do not match.")}
    }
  }

  ## Now check the dendrogram clusters and formatting, if dendrogram is provided
  if(!is.null(dend)){
    ## Check that dend is of class dendrogram
    if(!is.element("dendrogram",class(dend))){stop("If provided, dend must be of R class dendrogram.")}

    ## Check that dend and meta.data match in labels
    extra_labels <- setdiff(labels(dend), unique(as.character(meta.data[,celltypeColumn])))
    if(length(extra_labels)>0){stop(paste("Dendrogram has labels not included in metadata:",paste(extra_labels,collapse=", ")))}

    ## Check that meta.data and dend match in labels
    extra_labels <- setdiff(unique(as.character(meta.data[,celltypeColumn])), labels(dend))
    if(length(extra_labels)>0){
      warning(paste0("Metadata include cluster labels not found in dendrogram: ", paste(extra_labels, collapse=", ")))
    }
  }

  ## Check that cluster stats conforms to as.character(meta.data[,celltypeColumn])
  if(!is.null(cluster_stats)){
    print("===== Checking cluster_stats =====")
    if(!is.null(counts)){
      if(!all(rownames(cluster_stats) %in% colnames(counts))){
        stop("Cluster stats must have the same columns as the normalized expression matrix.")
      }
    }
    if(!all(unique(meta.data[[celltypeColumn]]) %in% colnames(cluster_stats))){
      stop("Cluster stats must have the same set of cell types as meta.data[[celltypeColumn]].")
    }
  }

  ## Ensure directory exists, if not create it
  taxonomyDir <- file.path(taxonomyDir) # Convert from windows to unix or vice versa
  dir.create(taxonomyDir, showWarnings = FALSE)

  ## Validate the schema
  # .schemaValidate(meta.data)

}

#' Function to subsample cells from a taxonomy
#'
#' @param cluster.names A vector of cell type names 
#' @param cell_ids A vector of cell ids to subsample from
#' @param dend A dendrogram object to use for subsampling
#' @param subsample The number of cells to retain per cluster (default = 2000)
#' 
#' @return boolean vector for subsampling
#'
#' @keywords internal
subsample_taxonomy = function(cluster.names, cell_ids, dend=NULL, subsample=2000){
  
  kpClusters <- rep(TRUE,length(cluster.names))
  if(!is.null(dend)){
    kpClusters <- is.element(cluster.names, labels(dend)) # exclude missing clusters, if any
    if(mean(kpClusters)<1) print("===== Omitting cells from missing clusters =====")
  }
  
  if((subsample > 0) & (subsample < Inf)){
      print("===== Subsampling cells =====")
      kpSub = cell_ids[subsampleCells(cluster.names, subsample)&kpClusters]
  }else{
      print("===== No subsampling of cells =====")
      kpSub = cell_ids[kpClusters]
  }

  return(kpSub)
}


#' Function to subsample cells
#'
#' @param cluster.names A vector of cluster names in the reference taxonomy.
#' @param subSamp Number of cells to keep per cluster.
#' @param seed Random seed used for subsampling.
#' @param use.historical Defualt (FALSE) uses new, faster implementation. Set to TRUE to use the historical, but slower implementation for back-compatibility.
#'
#' @return Boolean vector of cells to keep (TRUE) and cells to remove (FALSE)
#' 
#' @export
subsampleCells <- function(cluster.names, subSamp = 25, seed = 5, use.historical=FALSE) { 
  if(use.historical){
    kpSamp = rep(FALSE,length(cluster.names))
    for (cli in unique(as.character(cluster.names))){
      set.seed(seed)
      seed   = seed+1
      kp     = which(cluster.names==cli)
      kpSamp[kp[sample(1:length(kp),min(length(kp),subSamp))]] = TRUE
    }
    return(kpSamp)
  }
  # If use.historical is false, the new implementation (below) will be used
  if(length(subSamp)==1) 
    subSamp = rep(subSamp,length(unique(as.character(cluster.names)))) 
  if(is.null(names(subSamp))) 
    names(subSamp) <- unique(as.character(cluster.names)) 
  
  set.seed(seed) 
  cluster_split <- split(seq_along(cluster.names), as.character(cluster.names)) 
  
  kpSamp <- unlist(lapply(names(cluster_split), function(cli) { 
    val <- subSamp[cli] 
    if (!is.na(val)[1]) { 
      kp <- cluster_split[[cli]] 
      kp[sample(length(kp), min(length(kp), val))] 
    } else { 
      integer(0) 
    } 
  })) 
  kpSamp2 <- rep(FALSE, length(cluster.names)) 
  kpSamp2[kpSamp] <- TRUE
  kpSamp2
}


#' Function to set mapping mode
#'
#' @param AIT.anndata A vector of cluster names in the reference taxonomy.
#' @param mode Number of cells to keep per cluster.
#'
#' @return AIT anndata with mode set for mapping
#' 
#' @export
mappingMode <- function(AIT.anndata, mode){
  # Copied from scrattch.mapping
  if(!mode %in% names(AIT.anndata$uns$filter)){ stop(paste0(mode, " is invalid. Choose from: ", names(AIT.anndata$uns$filter))) }
  AIT.anndata$uns$mode = mode
  return(AIT.anndata)
}


#' Convert R dendrogram to json
#'
#' @param dend R dendrogram object
#'
#' @return json file
#' 
#' @export
dend_to_json = function(dend){
    # Convert dendrogram to hclust
    hclust_obj <- as.hclust(dend)
    ## Record information in list
    dendrogram_json <- list(
      node_heights = hclust_obj$height,
      cluster_tree = hclust_obj$merge,
      order = hclust_obj$order,
      labels = hclust_obj$labels)
    if(!all(is.na(get_nodes_attr(dend, "markers")))){
      dendrogram_json[["markers_names"]] = lapply(get_nodes_attr(dend, "markers"), names)
      dendrogram_json[["markers_values"]] = get_nodes_attr(dend, "markers")
      ## markers.byCl_names = lapply(get_nodes_attr(dend, "markesr.byCl"), names),
      ## markers.byCl_values = get_nodes_attr(dend, "markesr.byCl")
    }
    ## Return json
    return(dendrogram_json)
}


#' Convert json to R dendrogram
#'
#' @param json json from R dendrogram
#'
#' @return R dendrogram object
#' 
#' @export
json_to_dend = function(json){
    ## 
    hclust.tmp <- list()  # initialize empty object
    # define merging pattern: 
    #    negative numbers are leaves, 
    #    positive are merged clusters (defined by row number in $merge)
    json = fromJSON(json)
    hclust.tmp$merge <- json$cluster_tree    # leaf merges
    hclust.tmp$height <- json$node_heights   # define merge heights
    hclust.tmp$order <- json$order           # order of leaves(trivial if hand-entered)
    hclust.tmp$labels <- json$labels         # labels of leaves
    class(hclust.tmp) <- "hclust"                           # make it an hclust object
    dend = as.dendrogram(hclust.tmp)              # Make it an dendrogram object
    ##
    if("markers_names" %in% names(json)){
      ## Extract marker information per node from json data.frame
      markers = list()
      for(elem in 1:length(json$markers_values)){
          markers[[elem]] = setNames(json$markers_values[[elem]], json$markers_names[[elem]])
      }
      ## Add in marker information to dendrogram object
      i_node = 0
      get_attr_from_node <- function(dend_node, markers) {
          i_node <<- i_node + 1
          attr(dend_node, "markers") = markers[[i_node]]
          dend_node
      }
      ## Apply markers per node in a specific order based on dendextend
      dend = dendrapply(dend, get_attr_from_node, markers)
    }
    ##
    dend = label_dend(dend)[[1]]
    ## Return dendrogram
    return(dend)
}


#' Get top genes by beta (binary) score
#'
#' @param data A count (or CPM or logCPM) matrix
#' @param cluster.names A vector of cluster names in the reference taxonomy.
#' @param gene.count The number of top genes to return (Default=2000)
#'
#' @return Boolean vector of cells to keep (TRUE) and cells to remove (FALSE)
#'
#' @export
top_binary_genes <- function(data, cluster.names, gene.count=2000){
  cluster.names <- setNames(as.factor(cluster.names),colnames(data))
  propExpr  <- get_cl_prop(data,cluster.names)
  betaScore <- getBetaScore(propExpr,returnScore=FALSE)
  betaScore <- sort(betaScore)
  top.genes <- names(betaScore)[1:gene.count]
  return(top.genes)
}


####################################################################
## Functions for reversing '\' and '/'

#' Sets the default leading_string for file.path()
#' 
#' @param leading_string Default (NULL) sets to "\\\\" for Windows and "/" otherwise; or can provide any character vector
#'
#' @return leading_string for use as default in file.path().
#'
#' @export
setLeadingString <- function(leading_string=NULL){
  if (is.null(leading_string)){
    leading_string="/"
    if(get_os()=="windows") leading_string="\\\\" 
  }
  options("leading_string"=leading_string)
}


#' Sets the default path_separator for file.path()
#' 
#' @param path_separator Default (NULL) sets to .Platform$file.sep or can provide any character vector
#'
#' @return path_separator for use as default in file.path().
#'
#' @export
setPathSeparator <- function(path_separator=NULL){
  if (is.null(path_separator)){
    path_separator = .Platform$file.sep
  }
  options("path_separator"=path_separator)
}


#' Detect the operating system
#' 
#' This function was taken directly from https://conjugateprior.org/2015/06/identifying-the-os-from-r/ and all credit goes to Will Lowe from "conjugateprior".
#' 
#' @keywords internal
get_os <- function(){
  sysinf <- Sys.info()
  if (!is.null(sysinf)){
    os <- sysinf['sysname']
    if (os == 'Darwin')
      os <- "osx"
  } else { ## mystery machine
    os <- .Platform$OS.type
    if (grepl("^darwin", R.version$os))
      os <- "osx"
    if (grepl("linux-gnu", R.version$os))
      os <- "linux"
  }
  tolower(os)
}

#' Construct Path to File across platforms
#'
#' Construct the path to a file from components in a platform-independent way. This version is a wrapper of the base `file.path` function which also reverses '/' direction and can attempt to add double slashes if needed.
#'
#' @param ... character vectors. Long vectors are not supported.
#' @param path_separator the path separator to use (assumed to be ASCII).
#' @param leading_string what is the leading character(s) (e.g., '/' or '\\\\'). A string can be provided, or by default if the "leading_string" global variable is set it takes that variable, otherwise this is guessed at based on operating system and/or existence of a file at the file path
#' @param change_path_chars which characters should be changes to the path_separator value (default NULL = "none"). If you want to change all slashes in the file path to the correct direction set change_chars = c("/","\\"))
#' @param change_leading_chars which characters should be changes to the leading_string value (default = c("/","\\","\\\\")). This will preserve local file paths. 
#'
#' @return A file path with slashes going the correct direction.
#'
#' @export
file.path <- function (...,path_separator = getOption("path_separator"),leading_string=getOption("leading_string"),
                       change_path_chars=NULL,change_leading_chars=c("/","\\","\\\\")){
  if(is.null(path_separator))
    path_separator = .Platform$file.sep
  path <- base::file.path(...,path_separator)
  first_character <- grep('[^[:punct:]]', strsplit(path,"")[[1]])[1]
  start_path <- substring(path,1,(first_character-1))
  path <- substring(path,first_character,nchar(path))
  if (!is.null(change_path_chars)){
    path  <- strsplit(path,"")[[1]]
    slash <- which(is.element(path,change_path_chars)) 
    path[slash] <- path_separator
    path <- paste0(path,collapse="")
  }
  if (is.null(leading_string)){
    leading_string="/"
    if(get_os()=="windows") leading_string="\\\\" 
  }
  #remove trailing slashes
  while((substring(path,nchar(path),nchar(path)))==path_separator){
    path = substring(path,1,nchar(path)-1)
  }
  #add leading_string if desired
  add_lead = FALSE
  if(length(change_leading_chars)<1) 
    change_leading_chars="*"
  for (ch in change_leading_chars) if (grepl(ch,start_path, fixed=TRUE))
    add_lead = TRUE
  if(add_lead)
    path = paste0(leading_string,path)
  path
}

##################################################################################################################
## The functions below are mapping function from scrattch.hicat dev_zy branch that are required for tree mapping
##################################################################################################################

#' Function for building the standard reference format, including adding marker genes to the clustering tree
#'
#' @param cl Factor vector where values are cluster ids (e.g., a numeric vector of corresponding to cell type order in the tree) and values are sample ids for cells (e.g., this vector has length = number of cells) 
#' @param norm.dat log normalized expression data
#' @param dend Input dendrogram 
#' @param de.genes output from `display_cl` function
#' @param cl.label Factor vector here values are cluster ids (e.g., a numeric vector of corresponding to cell type order in the tree) and values are dendrogram labels (e.g., this vector has length = number of clusters) 
#' @param up.gene.score Output from `get_gene_score`
#' @param down.gene.score Output from `get_gene_score`
#' @param n.markers Number of marker genes to return per comparison (default=30)
#'
#' @return A list where `dend` is the updated dendrogram with markers attached and `cl.dat` is a matrix of cluster means
#'
#' @keywords internal
build_reference <- function(cl, norm.dat, dend, de.genes, cl.label=NULL, up.gene.score=NULL, down.gene.score=NULL, n.markers=30)
{
  suppressPackageStartupMessages({
    library(randomForest)
    library(scrattch.hicat)
  })
  
  cl.dat = get_cl_means(norm.dat, cl)
  if(is.null(up.gene.score)){
    de.gene.score = get_gene_score(de.genes)
    up.gene.score = de.gene.score[[1]]
    down.gene.score = de.gene.score[[2]]
  }    
  select.genes = intersect(row.names(norm.dat), row.names(up.gene.score))
  dend = select_dend_markers(dend, norm.dat=norm.dat, cl=cl, de.genes=de.genes,
                             up.gene.score=up.gene.score[select.genes,], 
                             down.gene.score=down.gene.score[select.genes,], n.markers=n.markers)
  dend = select_pos_dend_markers(dend= dend, norm.dat = norm.dat, cl = cl)
  if(!is.null(cl.label)){
    colnames(cl.dat) = cl.label[colnames(cl.dat)]
    labels(dend) = cl.label[labels(dend)]
  }
  dend = label_dend(dend)[[1]]
  labels(dend) <- setNames(colnames(cl.dat),colnames(cl.dat)) # This line might not be needed
  
  if(sum(!is.na(get_nodes_attr(dend, "original_label"))>0)){
    print("This section is needed if the starting dendrogram includes ccn nomenclature.")
    dend <- revert_dend_label(dend,get_nodes_attr(dend, "original_label"),"label")
  }
  return(list(cl.dat=cl.dat, dend=dend))
}

#' Compute cluster sums for each row in a matrix
#' 
#' This is the scrattch.hicat version of this function (the scrattch.bigcat version crashes the code).
#' 
#' @param mat A gene (rows) x samples (columns) sparse matrix
#' @param cl A cluster factor object
#' 
#' @return a matrix of genes (rows) x clusters (columns) with sums for each cluster
#' 
#' @keywords internal
get_cl_sums <- function(mat, 
                        cl) {
  
  cl.mat <- get_cl_mat(cl)
  if(all(names(cl) %in% colnames(mat))){
    cl.sums <- Matrix::tcrossprod(mat[,rownames(cl.mat)], Matrix::t(cl.mat))
  }
  else{
    cl.sums <- Matrix::crossprod(mat[rownames(cl.mat),], cl.mat)
  }
  cl.sums <- as.matrix(cl.sums)
  return(cl.sums)
}

#' Compute cluster means for each row in a matrix
#' 
#' This is the scrattch.hicat version of this function (the scrattch.bigcat version crashes the code).
#' 
#' @param mat A gene (rows) x samples (columns) sparse matrix
#' @param cl A cluster factor object
#' 
#' @return a matrix of genes (rows) x clusters (columns) with means for each cluster
#' 
#' @keywords internal
get_cl_means <- function (mat, cl) 
{
  cl.sums <- get_cl_sums(mat, cl)
  cl.size <- table(cl)
  cl.means <- as.matrix(Matrix::t(Matrix::t(cl.sums)/as.vector(cl.size[colnames(cl.sums)])))
  return(cl.means)
}

#' Strip extra annotation information from dendrogram
#'
#' @param dend R dendrogram object
#' @param value Vector of values pulled from the dendrogram
#' @param attribute Which attribute should be overwritten
#'
#' @return R dendrogram object with updated attributes
#'
#' @keywords internal
revert_dend_label <- function(dend, value, attribute="label")
{
  if(attr(dend, attribute)=="")
    attr(dend, attribute) <- value[attr(dend,"original_label")]
  if (length(dend)>1) for(i in 1:length(dend))
    dend[[i]]=revert_dend_label(dend[[i]], value=value, attribute)
  return(dend)
}

#' map_dend_membership
#'
#' @param dend R dendrogram in a specific format
#' @param cl.dat gene by cell type matrix (I think?)
#' @param map.dat normalized data of the MAPPING data set.
#' @param map.cells names of cells to map (e.g., the column names of the cell x gene matrix)
#' @param mc.cores number of cores to run the mapping on 
#' @param bs.num Number of bootstrapping runs to calculate the membership from (default = 100)
#' @param seed = random seed
#' @param ... other variables to pass to map_dend
#'
#' @import foreach
#'
#' @return membership table
#' 
#' @keywords internal
.map_dend_membership <-
  function(dend,
           cl.dat,
           map.dat,
           map.cells,
           mc.cores = 10,
           bs.num = 100,
           seed = 42,
           ...)
  {
    
    # Optional libraries for UNIX parallel implementation (likely will crash in Windows)
    if (mc.cores>1) {
      suppressPackageStartupMessages({
        library(doMC)
        library(parallel)
      })
    } 
    
    if(mc.cores ==1){
      registerDoSEQ()
    }else{
      registerDoMC(cores=mc.cores)
    }
    mem = foreach(i = 1:bs.num, .combine = 'c') %dopar% {
      print(i)
      .map_dend(dend, cl.dat, map.dat, map.cells, seed=i)
    }
    memb = data.frame(cell = names(mem), cl = mem)
    memb = table(memb$cell, memb$cl)
    memb = memb / bs.num
    tmp = get_nodes_attr(dend, "label")
    tmp = tmp[tmp %in% colnames(memb)]
    memb = memb[, tmp]
    return(memb)
  }

#' map_dend
#'
#' @param dend A dendrogram in R format 
#' @param cl A cluster factor object to compare to a reference
#' @param cl.dat gene by cell type matrix (I think?)
#' @param map.dat normalized data of the MAPPING data set.
#' @param select.cells names of cells to map (e.g., the column names of the cell x gene matrix)
#' @param p proportion of marker genes to include in each iteration of the mapping algorithm.
#' @param low.th the minimum difference in Pearson correlation required to decide on which branch
#' @param default.markers What genes to include in every bootstrap run (default is none)
#' @param seed = random seed
#' 
#' @return tree mapping to the dendrogram table (cells x nodes with values as probabilities)
#' 
#' @keywords internal
.map_dend <-
  function(dend,
           cl.dat,
           map.dat,
           select.cells=colnames(map.dat),
           p = 0.8,
           low.th = 0.1,
           default.markers = NULL,
           seed = 42)
  {
    final.cl = c(setNames(rep(
      attr(dend, "label"), length(select.cells)
    ), select.cells))
    if (length(dend) <= 1) {
      return(final.cl)
    }
    markers = attr(dend, "markers")
    markers = markers[names(markers) %in% row.names(map.dat)]
    cl.g = sapply(dend, labels, simplify = F)
    names(cl.g) = 1:length(cl.g)
    genes = names(markers)
    genes = union(genes, default.markers)
    mapped.cl = .resolve_cl(cl.g,
                           cl.dat,
                           markers,
                           map.dat,
                           select.cells,
                           p = p,
                           low.th = low.th,
                           seed = seed+1)
    if (length(mapped.cl) > 0) {
      for (i in unique(mapped.cl)) {
        select.cells = names(mapped.cl)[mapped.cl == i]
        if (length(select.cells) > 0) {
          final.cl = c(
            final.cl,
            .map_dend(
              dend[[as.integer(i)]],
              cl.dat,
              map.dat,
              select.cells,
              p = p,
              low.th = low.th,
              seed = seed+2
            )
          )
        }
      }
      return(cl = final.cl)
    }
    
  }

#' resolve_cl
#'
#' @param cl.g Cluster labels in some format
#' @param cl.med Cluster medians
#' @param markers Genes to use as markers for this function
#' @param map.dat normalized data of the MAPPING data set.
#' @param select.cells names of cells to map (e.g., the column names of the cell x gene matrix)
#' @param p proportion of marker genes to include in each iteration of the mapping algorithm.
#' @param low.th the minimum difference in Pearson correlation required to decide on which branch
#' @param seed - random seed for reproducibility
#'
#' @return mapped.cl output
#' 
#' @keywords internal
.resolve_cl <-
  function(cl.g,
           cl.dat,
           markers,
           map.dat,
           select.cells,
           p = 0.8,
           low.th = 0.1,
           seed = 42)
  {
    ##
    genes = names(markers)[markers > 0]
    tmp.cl = unlist(cl.g)
    
    ###For each branch point, find the highest expression cluster.
    tmp.med = sapply(cl.g, function(g) rowMaxs(cl.dat[genes, g, drop = F]))
    row.names(tmp.med) = genes
    ###Make sure the genes are discriminative between all the branches.
    genes = genes[rowMaxs(tmp.med) - rowMins(tmp.med) > 1]
    
    ###Sample the markers based on the weights.
    ##TO DO: randomforest sometimes give importance value of 0. adjust for that.
    set.seed(seed)
    seed  = seed+1
    genes = sample(genes, round(length(genes) * p), prob = markers[genes])
    
    ###Compute the correlation with the median cluster profile.
    ###add drop=F
    cl.cor = cor(as.matrix(map.dat[genes, select.cells, drop = F]), cl.dat[genes, tmp.cl, drop =
                                                                             F])
    cl.cor[is.na(cl.cor)] = 0
    ###Compute the best match in each branch.
    tmp.score = do.call("cbind", sapply(cl.g, function(x)
      rowMaxs(cl.cor[, x, drop = F]), simplify = F))
    row.names(tmp.score) = row.names(cl.cor)
    ####Determine the best match.
    best.score = setNames(rowMaxs(tmp.score), row.names(tmp.score))
    ###determine the difference from the best match.
    diff.score = best.score - tmp.score
    
    ####Give up on cells can't be discriminated,choose one branch randomly.
    unresolved.cl = row.names(tmp.score)[rowSums(diff.score < low.th) ==
                                           ncol(diff.score)]
    set.seed(seed)
    seed  = seed+1
    mapped.cl = setNames(sample(colnames(tmp.score), length(unresolved.cl), replace =
                                  T), unresolved.cl)
    
    ###Cells mapped to one or more branches.
    mapped.cells = setdiff(row.names(cl.cor), unresolved.cl)
    ###For binary branch, done already
    if (length(cl.g) == 2) {
      mapped.cl = c(mapped.cl, setNames(colnames(diff.score)[apply(diff.score[mapped.cells, , drop =
                                                                                F], 1, which.min)], mapped.cells))
      return(mapped.cl)
    }
    ##The remaining options for mapped cells
    tmp.cl = sapply(mapped.cells, function(x)
      colnames(diff.score)[which(diff.score[x,] < low.th)], simplify = F)
    ###cells with multiple options
    resolve.cells = names(tmp.cl)[sapply(tmp.cl, length) > 1]
    ###cells with only one option. Not further job.
    mapped.cells = setdiff(mapped.cells, resolve.cells)
    if (length(mapped.cells) > 0) {
      mapped.cl = c(mapped.cl, setNames(unlist(tmp.cl[mapped.cells]), mapped.cells))
    }
    ###Resolve further options.
    if (length(resolve.cells) > 0) {
      tmp.cat = sapply(tmp.cl[resolve.cells], function(x)
        paste(x, collapse = " "))
      for (cat in unique(tmp.cat)) {
        tmp.cl = unlist(strsplit(cat, " "))
        select.cells = names(tmp.cat)[tmp.cat == cat]
        mapped.cl = c(
          mapped.cl,
          resolve_cl(
            cl.g[tmp.cl],
            cl.dat,
            markers,
            map.dat,
            select.cells,
            p = p,
            low.th = low.th
          )
        )
      }
    }
    return(mapped.cl)
  }

#' Build dend (updated to specify dendextend version of "set")
#'
#' @param cl.dat Normalized data of the REFERENCE data set
#' @param cl.cor Matrix of cell x cell correlations (calculated if not provided) 
#' @param l.rank Factor of cluster order (in a specific format)
#' @param l.color Factor of clluster colors (in a specific format)
#' @param nboot Number of bootstrapping runs to calculate the membership from (default = 100)
#' @param ncores Number of cores for performing calculations
#'
#' @return dendrogram and a couple of related things
#'
#' @import dendextend
#' @import pvclust
#' 
#' @keywords internal
build_dend <- function(cl.dat, cl.cor=NULL, l.rank=NULL, l.color=NULL, nboot=100, ncores=1)
{
  if(is.null(cl.cor)){
    cl.cor = cor(cl.dat)
  }
  pvclust.result=NULL
  if(nboot > 0){
    require(pvclust)
    parallel= FALSE
    if(ncores > 1){
      parallel = as.integer(ncores)
    }
    pvclust.result <- pvclust::pvclust(cl.dat, method.dist = "cor" ,method.hclust = "average", nboot=nboot, parallel=parallel)
    dend = as.dendrogram(pvclust.result$hclust)
    dend = label_dend(dend)$dend
    dend = dend %>% pvclust_show_signif_gradient(pvclust.result, signif_type = "bp", signif_col_fun=colorRampPalette(c("white","gray","darkred","black")))
  }
  else{
    cl.hc = hclust(as.dist(1-cl.cor),method="average")      
    dend = as.dendrogram(cl.hc)
  }
  dend = dend %>% dendextend::set("labels_cex", 0.7)
  if(!is.null(l.color)){
    dend = dend %>% dendextend::set("labels_col", l.color[labels(dend)])
  }
  dend = dend %>% dendextend::set("leaves_pch", 19) %>% dendextend::set("leaves_cex", 0.5)
  if(!is.null(l.color)){
    dend = dend %>% dendextend::set("leaves_col", l.color[labels(dend)])
  }
  if(!is.null(l.rank)){
    dend =reorder_dend(dend,l.rank)
  }
  return(list(dend=dend, cl.cor=cl.cor, pvclust.result=pvclust.result))
}

#' Compute cluster medians for each row in a matrix
#' 
#' @param mat A gene (rows) x samples (columns) sparse matrix
#' @param cl A cluster factor object
#' 
#' @return a matrix of genes (rows) x clusters (columns) with medians for each cluster
#'
#' @export
get_cl_medians <- function(mat, cl)
{
  library(Matrix)
  library(matrixStats)
  
  cl.med <- do.call("cbind",
                    tapply(names(cl), 
                           cl, 
                           function(x){
                             matrixStats::rowMedians(as.matrix(mat[,x]))
                           }
                    )
  )
  
  rownames(cl.med) <- rownames(mat)
  
  return(cl.med)
}


#' Convert a matrix of raw counts to a matrix of Counts per Million values
#' 
#' The input can be a base R matrix or a sparse matrix from the Matrix package.
#' 
#' This function expects that columns correspond to samples, and rows to genes, but can also take the transpose if specified by cells.as.rows=TRUE.
#' 
#' @param counts a matrix, dgCMatrix, or dgTMatrix of count values.
#' @param sf vector of numeric values representing the total number of reads. If count matrix includes all genes, value calulated by default (sf=NULL) will be accurate; however, if count matrix represents only a small fraction of genes, we recommend also providing this value.
#' @param denom Denominator that all counts will be scaled to. The default (1 million) is commonly used, but 10000 is also common for sparser droplet-based sequencing methods.
#' @param cells.as.rows Set to FALSE (default) if rows are genes and columns are cells or TRUE if rows are cells and columns are genes.
#' 
#' 
#' @return a matrix, dgCMatrix, or dgTMatrix of CPM values (matching input)
#' 
#' @export 
cpm <- function (counts, 
                 sf = NULL, 
                 denom = 1e+06, 
                 cells.as.rows=FALSE) 
{
  if(cells.as.rows){
    # Cells as rows, genes as columns
    if((intersect(class(counts),c("dgTMatrix","dgCMatrix"))>0)|(is.matrix(counts))) {
      # Calculate total counts per cell (now per row)
      if (is.null(sf)) {
        sf <- Matrix::rowSums(counts)
      }
      sf = sf/denom
      
      # Calculate CPM for each gene in each cell
      cpm = sweep(counts, 1, sf, "/", check.margin = FALSE)
      
      ## FUTURE UPDATE: make more memory efficient for dgCMatrix and dgTMatrix classes like below
      
      return(cpm)
    } 
    else {
      stop(paste("cpm function for", class(counts)[1], "not supported"))
    }
  } 
  else { 
    # Cells as columns, genes as rows
    if (is.null(sf)) {
      sf <- Matrix::colSums(counts)
    }
    sf = sf/denom
    if (is.matrix(counts)) {
      return(sweep(counts, 2, sf, "/", check.margin = FALSE))
    }
    else if (class(counts) == "dgCMatrix") {
      sep <- counts@p
      sep <- sep[-1] - sep[-length(sep)]
      j <- S4Vectors::Rle(1:length(sep), sep)
      counts@x <- counts@x/sf[as.integer(j)]
    }
    else if (class(counts) == "dgTMatrix") {
      j = counts@j
      counts@x = counts@x/sf[j + 1]
    }
    else {
      stop(paste("cpm function for", class(counts)[1], "not supported"))
    }
    return(counts)
  }
}


#' Convert a matrix of raw counts to a matrix of log2(Counts per Million + 1) values
#' 
#' The input can be a base R matrix or a sparse matrix from the Matrix package.
#' 
#' This function expects that columns correspond to samples, and rows to genes by default, but can also take the transpose if specified by cells.as.rows=TRUE.  By default the offset is 1, but to calculate just log2(counts per Million) set offset to 0.
#' 
#' @param counts A matrix, dgCMatrix, or dgTMatrix of count values
#' @param offset The constant offset to add to each cpm value prior to taking log2 (default = 1)
#' @param sf vector of numeric values representing the total number of reads. If count matrix includes all genes, value calulated by default (sf=NULL) will be accurate; however, if count matrix represents only a small fraction of genes, we recommend also providing this value.
#' @param denom Denominator that all counts will be scaled to. The default (1 million) is commonly used, but 10000 is also common for sparser droplet-based sequencing methods.
#' @param cells.as.rows Set to FALSE (default) if rows are genes and columns are cells or TRUE if rows are cells and columns are genes.
#' 
#' @return a matrix, dgCMatrix, or dgTMatrix of log2(CPM + 1) values (matching input)
#' 
#' @export 
logCPM <- function (counts, 
                    offset = 1, 
                    sf = NULL, 
                    denom = 1e+06, 
                    cells.as.rows=FALSE,
                    ...) 
{
  # Run log2CPM_byRow if cells are rows and counts is sparse
  if((as.character(class(counts))=="dgCMatrix")&(cells.as.rows))
    return(log2CPM_byRow(counts=counts, sf=sf, denom=denom, offset=offset))
  
  # Otherwise run standard function
  norm.dat <- cpm(counts, sf=sf, denom=denom, cells.as.rows=cells.as.rows)
  if (is.matrix(norm.dat)) {
    norm.dat <- log2(norm.dat + offset)
  }
  else {
    norm.dat@x <- log2(norm.dat@x + offset)
  }
  norm.dat
}


#' Convert a matrix of raw counts to a matrix of log2(Counts per Million + 1) values
#' 
#' The input can be a base R matrix or a sparse matrix from the Matrix package.
#' 
#' This function expects that columns correspond to genes, and rows to samples by default and is equivalent to running logCPM with cells.as.rows=TRUE (but a bit faster).  By default the offset is 1, but to calculate just log2(counts per Million) set offset to 0.
#' 
#' @param counts A matrix, dgCMatrix, or dgTMatrix of count values
#' @param sf vector of numeric values representing the total number of reads. If count matrix includes all genes, value calulated by default (sf=NULL) will be accurate; however, if count matrix represents only a small fraction of genes, we recommend also providing this value.
#' @param denom Denominator that all counts will be scaled to. The default (1 million) is commonly used, but 10000 is also common for sparser droplet-based sequencing methods.
#' @param offset The constant offset to add to each cpm value prior to taking log2 (default = 1)
#' 
#' @return a dgCMatrix of log2(CPM + 1) values
#' 
#' @export 
log2CPM_byRow <- function (counts, sf = NULL, denom = 1e+06, offset=1){
  if(!(as.character(class(counts))=="dgCMatrix"))
    counts <- as(counts, "dgCMatrix")
  if (is.null(sf)) {
    sf <- Matrix::rowSums(counts)
  }
  sf <- sf/denom
  normalized.expr   <- counts
  normalized.expr@x <- counts@x/sf[as.integer(counts@i + offset)]
  normalized.expr@x <- log2(normalized.expr@x+1)
  return(normalized.expr)
}

