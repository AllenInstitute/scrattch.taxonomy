#' Function to define AIT schema
#'
#' @return AIT schema
#'
#' @keywords internal
.schemaAIT = function(){
  data = data.frame("")
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
#' @param feature.set A list of variable features
#' @param umap.coords A matrix of UMAP coordinates
#' 
#' @return Stops the function if any of the parameters are not as expected
#'
#' @keywords internal
.checkBuildTaxonomyParams <- function(counts, meta.data, feature.set, 
                                        umap.coords, taxonomyDir, taxonomyName, 
                                        celltypeColumn, cluster_colors, cluster_stats, dend){
  if(sum(is.element(paste0(celltypeColumn,c("","_label")), colnames(meta.data)))==0){stop("cluster column must be defined in the meta.data object")}
  if(is.null(feature.set)){stop("Compute variable features and supply feature.set")}
  if(is.null(umap.coords)){stop("Compute UMAP dimensions and supply umap.coords")}
  if(!all(colnames(counts) == rownames(meta.data))){stop("Colnames of `counts` and rownames of `meta.data` do not match.")}
  if(!is.data.frame(meta.data)){stop("meta.data must be a data.frame, convert using as.data.frame(meta.data)")}
  if("sample" %in% colnames(meta.data)){stop("meta.data column name 'sample' is reserved and cannot be used, please remove or rename.")}

  ## Rename celltypeColumn to "cluster" if needed
  celltypeColumn <- gsub("_label","",celltypeColumn)
  if(celltypeColumn!="cluster"){
    inCol     <- paste0(celltypeColumn,c("","_label","_id","_color"))
    outCol    <- paste0("cluster",c("","_label","_id","_color"))
    meta.data <- meta.data[,setdiff(colnames(meta.data),outCol)]
    for (i in 1:4) if(is.element(inCol[i],colnames(meta.data))) {
      colnames(meta.data) <- gsub(inCol[i],outCol[i],colnames(meta.data))
    }
  }
  if(sum(is.element(c("cluster","cluster_label"),colnames(meta.data)))>1){stop("Only a single cluster column can be provided (e.g., cluster or cluster_label but not both).")}
  
  ## Capture the cluster colors from the metadata if provided and if possible
  if(sum(is.element("cluster_color", colnames(meta.data))) == 1){
    if(length(meta.data$cluster_label)>0){
      cluster_colors <- setNames(meta.data$cluster_color, meta.data$cluster_label)
      cluster_colors <- cluster_colors[unique(names(cluster_colors))]
      meta.data$cluster <- meta.data$cluster_label
      meta.data <- meta.data[,setdiff(colnames(meta.data),paste0("cluster",c("_label","_id","_color")))]
    }else {
      warning("Cannot match cluster_label and cluster_color in meta.data, so cluster_color will be ignored.")
    }
  }

  ## Now check the dendrogram clusters and formatting, if dendrogram is provided
  if(!is.null(dend)){
    if(!is.element("dendrogram",class(dend))){stop("If provided, dend must be of R class dendrogram.")}
    clusters=unique(meta.data$cluster)
    extra_labels <- setdiff(labels(dend), clusters)
    if(length(extra_labels)>0){stop(paste("Dendrogram has labels not included in metadata:",paste(extra_labels,collapse=", ")))}
    extra_labels <- setdiff(clusters, labels(dend))
    if(length(extra_labels)>0){
      warning(paste0("Metadata include cluster labels not found in dendrogram: ", paste(extra_labels, collapse=", "),
                     ". Cells from these clusters will be EXCLUDED from all taxonomy files."))
    }
  }

  ## Check that cluster stats conforms to meta.data$cluster
  if(!is.null(cluster_stats)){
    if(!all(unique(meta.data$cluster) %in% cluster_stats$cluster)){
      stop("clusters in provided cluster_stats must contain cover all of meta.data$cluster")
    }
  }

  ## Ensure directory exists, if not create it
  taxonomyDir <- file.path(taxonomyDir) # Convert from windows to unix or vice versa
  dir.create(taxonomyDir, showWarnings = FALSE)

  ##
  return(list(meta.data=meta.data, celltypeColumn=celltypeColumn))
}

#' Function to subsample cells
#'
#' @param cluster.names A vector of cluster names in the reference taxonomy.
#' @param subSamp Number of cells to keep per cluster.
#' @param seed Random seed used for subsampling.
#'
#' @return Boolean vector of cells to keep (TRUE) and cells to remove (FALSE)
#' 
#' @export
subsampleCells <- function(cluster.names, subSamp=25, seed=5){
  # Returns a vector of TRUE false for choosing a maximum of subsamp cells in each cluster
  # cluster.names = vector of cluster labels in factor format
  kpSamp = rep(FALSE,length(cluster.names))
  for (cli in unique(as.character(cluster.names))){
    set.seed(seed)
    seed   = seed+1
    kp     = which(cluster.names==cli)
    kpSamp[kp[sample(1:length(kp),min(length(kp),subSamp))]] = TRUE
  }
  return(kpSamp)
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
#' @keywords external
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