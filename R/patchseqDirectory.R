#' Starting from an anndata object this function builds the minimum files required for patch-seq shiny
#'
#' @param AIT.anndata A reference taxonomy object.
#' @param mappingFolder The location to save output files for patch-seq (or other query data) results, e.g. "/allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/star/human/human_patchseq_MTG_JAM_TEST/current/"
#' @param query.data A CPM normalized matrix to be annotated.
#' @param query.metadata A data frame of metadata for the query data.  
#' @param query.mapping Mapping results from `taxonomy_mapping()`, must be an mappingClass S4 object. If provided row names must match column names in query.data.
#' @param doPatchseqQC Boolean indicating whether patch-seq QC metrics should be calculated (default) or not.
#' @param metadata_names An optional named character vector where the vector NAMES correspond to columns in the metadata matrix and the vector VALUES correspond to how these metadata should be displayed in Shiny. This is used for writing the desc.feather file later.
#' @param min.confidence Probability below which a cell cannot be assigned to a cell type (default 0.7).  In other words, if no cell types have probabilities greater than resolution.index, then the assigned cluster will be an internal node of the dendrogram. 
#' @param verbose Should detail logging be printed to the screen?
#' 
#' This function writes files to the mappingFolder directory for visualization with molgen-shiny tools  
#' --- anno.feather - query metadata  
#' --- data.feather - query data  
#' --- dend.RData   - dendrogram (copied from reference)  
#' --- desc.feather - table indicating which anno columns to share  
#' --- memb.feather - tree mapping of each query cell to each tree node (not just the best matching type like in treeMap)  
#' --- tsne.feather - low dimensional coordinates for data  
#' --- tsne_desc.feather - table indicating which low-D representations to share  
#' 
#' @import reticulate
#' 
#' @export
buildMappingDirectory = function(AIT.anndata, 
                                 mappingFolder,
                                 query.data,
                                 query.metadata,
                                 query.mapping,
                                 doPatchseqQC = TRUE,
                                 metadata_names = NULL,
                                 min.confidence = 0.7,
                                 verbose=TRUE
){

  ## Gather mapping results from S4 class
  mapping.results = getMappingResults(query.mapping)
  
  ## Checks
  if(!all(colnames(query.data) == rownames(query.metadata))){stop("Colnames of `query.data` and rownames of `query.metadata` do not match.")}
  if(!all(colnames(query.data) == rownames(mapping.results))){stop("Colnames of `query.data` and rownames of `mapping.results` do not match.")}

  ## Filter taxonomy to mode cells
  AIT.anndata = AIT.anndata[!AIT.anndata$uns$filter[[AIT.anndata$uns$mode]]]  # This is a filter of cells to OMIT

  ## Merge mapping metadata inputs
  if(length(setdiff(colnames(mapping.results),colnames(query.metadata)))>0){
    query.metadata <- cbind(query.metadata, mapping.results[,setdiff(colnames(mapping.results), colnames(query.metadata))])
    query.metadata[,colnames(mapping.results)] <- mapping.results # Overwrite any columns in query.metadata with same name in mapping.results
  }
  
  ## Ensure directory exists, if not create it
  mappingFolder <- file.path(mappingFolder) ## Allow for unix or windows
  dir.create(mappingFolder, recursive=T, showWarnings = FALSE)

  ##
  if(!file.exists(mappingFolder)) {stop(paste("Directory ",mappingFolder," could not be created. Please do so manually or try again."))}

  if(verbose == TRUE) print("Gathering cluster medians from taxonomy folder.")

  ## Read in cluster medians from uns file location (either within anndata object or from the folder [for back-compatibility])
  if(!is.null(AIT.anndata$uns$stats[["standard"]]$medianmat)){
    cl.summary = reticulate::py_to_r(AIT.anndata$uns$stats[["standard"]]$medianmat)  # Medians now stored in the uns to avoid needing to read files
  } else {
    cl.summary = read_feather(file.path(AIT.anndata$uns$taxonomyDir, "medians.feather")) %>% as.data.frame()
  }
  cl.dat = as.matrix(cl.summary[,-1]); rownames(cl.dat) = cl.summary[,1]

  if(verbose == TRUE) print("Saving dendrogram to mapping folder.")

  ## Read in the reference tree and copy to new directory
  dend = json_to_dend(AIT.anndata$uns$dend[[AIT.anndata$uns$mode]])

  ## Output dend to mapping folder
  saveRDS(dend, file.path(mappingFolder, "dend.RData"))
  
  ## Convert query.data to CPM
  if(is.element("data.frame",class(query.data))){stop("`query.data` should be a matrix or a sparse matrix, not a data.frame.")}
  query.cpm <- cpm(query.data)
  sample_id <- colnames(query.cpm); query.metadata$sample_id = sample_id
  gene      <- rownames(query.cpm)

  ##
  binary.genes <- intersect(AIT.anndata$var_names[AIT.anndata$var$highly_variable_genes], rownames(query.cpm))

  ## Check for cells with empty data
  bad.cells    <- which(colSums(query.cpm[binary.genes,]>0)<=1)
  if(length(bad.cells)>0){
    query.cpm[,bad.cells] <- rowMeans(query.cpm)
    warning(paste("WARNING: the following query cells do not express any marker genes and are almost definitely bad cells:",
                  paste(colnames(query.cpm)[bad.cells],collapse=", ")))
  }

  if(verbose == TRUE) print("Saving tree mapping membership table to mapping folder.")
  
  ## Create and output the memb.feather information
  memb.ref <- query.mapping@detailed_results[["tree"]]
  memb.ref <- memb.ref[sample_id,]
  memb     <- data.frame(sample_id,as.data.frame.matrix(memb.ref)) 
  colnames(memb) <- c("sample_id",colnames(memb.ref))
  write_feather(as_tibble(memb), file.path(mappingFolder, "memb.feather")) 
  
  ## Assign lowest node (or leaf) on tree with confidence > min.confidence to the metadata "cluster" column
  confident_match <- apply(memb,1,function(x){
    i = max(which(x>=min.confidence))
    c(colnames(memb)[i],x[i])
  })
  confident_match <- as.data.frame(t(confident_match))
  ## getTopMatch(memb.ref[,intersect(labels(dend),colnames(memb.ref))]) # If we want to require leaf node, this line replaces above
  colnames(confident_match) <- c("cluster", "cluster_score")
  query.metadata <- cbind(query.metadata[,setdiff(colnames(query.metadata),c("cluster","cluster_score"))],confident_match[sample_id,])
  
  ## Define resolution index and cluster order, and update metadata
  invisible(capture.output({  # Avoid printing lots of numbers to the screen
    ## Get shiny's cluster_ids for each node
    cl.order = dend %>% get_nodes_attr("label")
    cluster_anno_ordered = data.frame(cluster_label=cl.order,
                                      cluster_id=1:length(cl.order))
    rownames(cluster_anno_ordered) = cl.order
    nodes.ids.df = data.frame(cluster_height = get_nodes_attr(dend,"height"),
                              cluster_label = get_nodes_attr(dend,"label"))
    nodes.ids.df$res.index = 1 - (nodes.ids.df$cluster_height / attr(dend,"height"))
    ## Join the memb with the cluster names and cluster res index
    cluster_node_anno <- cluster_anno_ordered %>%
      left_join(nodes.ids.df) %>%
      select(cluster_label, cluster_id, res.index)
    ## Ordered by dend and higher res first
    idx1 = which(cluster_node_anno$res.index>=0.8)
    idx2 = which(cluster_node_anno$res.index<0.8)
    idx3 = which(is.na(cluster_node_anno$res.index))
    cluster_node_anno = cluster_node_anno[c(idx1,idx2,idx3),]
    cluster_node_anno$cluster_id = 1:length(cluster_node_anno$cluster_label)
  }))

  ## Add this stuff to the meta.data file
  query.metadata$res.index <- cluster_node_anno$res.index[match(query.metadata$cluster,cluster_node_anno$cluster_label)]
  query.metadata$cluster   <- droplevels(factor(query.metadata$cluster, levels=cluster_node_anno$cluster_label))
  
  ## Output query data feather
  if(verbose == TRUE) print("Writing query data feather to mapping folder.")
  norm.data.t = Matrix::t(as.matrix(query.cpm))
  norm.data.t = as_tibble(norm.data.t)
  norm.data.t = cbind(sample_id, norm.data.t)
  norm.data.t = as_tibble(norm.data.t)
  write_feather(norm.data.t, file.path(mappingFolder, "data.feather"))
  
  ## Apply PatchseqQC if desired
  if(doPatchseqQC == TRUE){
    if(verbose == TRUE) print("Performing patchseqQC and updating metadata.")
    query.metadata <- applyPatchseqQC(AIT.anndata, query.data, query.metadata)
  }
  
  ## Add quality calls from KL divergence to the query metadata [***NEW!***]
  query.metadata <- cbind(query.metadata, tree_quality_call(AIT.anndata, query.mapping))
  
  ## Process and output metadata and desc files
  # Auto_annotate the data
  meta.data = scrattch.io::auto_annotate(query.metadata)
  
  # Convert chars and factors to characters (can potentially omit)
  for (col in colnames(meta.data)){ 
    if(is.character(meta.data[,col]) | is.factor(meta.data[,col])){
      meta.data[,col] = as.character(meta.data[,col])
    }
  }
  # Fix varibow color set
  for (col in which(grepl("_color",colnames(meta.data)))){
    kp = nchar(meta.data[,col])==5
    meta.data[kp,col] = paste0(meta.data[kp,col],"FF")
  }
  
  # Adjust the cluster colors to match cluster_colors in the reference taxonomy
  ln  <- dim(meta.data)[2]
  label.cols <- Reduce(intersect, list(c("cluster_label","subclass_label", "class_label"), colnames(meta.data), colnames(AIT.anndata$obs))) ## Check intersect against both meta.data and taxonomy. Edge case: query contains "class_label" and taxonomy does not.
  for(lab in label.cols){
    cnt <- setNames(rep(0,ln),colnames(meta.data))
    for (i in 1:ln) if(mean(is.element(meta.data[,i], (AIT.anndata$obs[,lab])))==1){
      col <- gsub("_label","",colnames(meta.data)[i])
      print(paste("Colors updated for:",col))
      meta.data[,paste0(col,"_color")] <- AIT.anndata$obs[,gsub("_label","_color",lab)][match(meta.data[,paste0(col,"_label")],AIT.anndata$obs[,lab])]
    }
  }
 
  ## Write the desc file.  
  anno_desc = create_desc(meta.data, use_label_columns = TRUE)
  # Subset the desc file to match metadata_names, if provided
  if(!is.null(metadata_names)){
    desc <- anno_desc[match(names(metadata_names), as.character(as.matrix(anno_desc[,1]))),]
    desc[,2] <- as.character(metadata_names[as.character(as.matrix(desc[,1]))])
    desc <- desc[!is.na(desc$base),]  # Remove missing values
    anno_desc <- desc
  }
  write_feather(as_tibble(anno_desc), file.path(mappingFolder, "desc.feather"))
  
  ## Minor reformatting of metadata file, then write metadata file
  if(verbose == TRUE) print("Writing annotation to mapping folder.")
  meta.data$cluster = meta.data$cluster_label; # Not sure why this is needed
  colnames(meta.data)[colnames(meta.data)=="sample_name"] <- "sample_name_old" # Rename "sample_name" to avoid shiny crashing
  if(!is.element("sample_id", colnames(meta.data))){ meta.data$sample_id = meta.data$sample_name_old } ## Sanity check for sample_id
  write_feather(meta.data, file.path(mappingFolder,"anno.feather"))
  
  ## Project mapped data into existing umap (if it exists) or generate new umap otherwise
  ref.umap     <- as.matrix(AIT.anndata$obsm[["umap"]][,colnames(AIT.anndata$obsm[["umap"]])!="sample_id"])
  rownames(ref.umap) <- rownames(AIT.anndata$obsm[["umap"]])
  ref.umap[is.na(ref.umap)] <- 0
  
  ##
  if(verbose == TRUE) print("Building UMAP for query data.")
  npcs         <- min(30,length(binary.genes))
  query.pcs    <- prcomp(logCPM(query.cpm)[binary.genes,], scale = TRUE)$rotation
  
  if(diff(range(ref.umap))>0){
    reference.logcpm <- t(AIT.anndata$X[,binary.genes])
    
    ## Check for cells with empty data
    bad.cells    <- which(colSums(reference.logcpm[binary.genes,]>0)<=1)
    if(length(bad.cells)>0){
      reference.logcpm[,bad.cells] <- rowMeans(reference.logcpm)
      warning(paste("WARNING: the following reference cells do not express any marker genes and are almost definitely bad cells:",
                    paste(colnames(reference.logcpm)[bad.cells],collapse=", ")))
    }
    
    ## Project mapped data into existing umap space
    reference.pcs    <- prcomp(reference.logcpm, scale = TRUE)$rotation
    reference.umap   <- umap(reference.pcs[,1:npcs])
    reference.umap$layout <- ref.umap[rownames(reference.umap$layout),]
    query.umap <- predict(reference.umap, query.pcs[,1:npcs])
  } else {
    ## Calculate a UMAP based on PCS of variable genes
    query.umap <- umap(query.pcs[,1:npcs])$layout
  }

  ## Write the UMAP coordinates.  
  tsne      = data.frame(sample_id = rownames(query.umap),
                         all_x = query.umap[,1],
                         all_y = query.umap[,2])
  tsne      = tsne[match(meta.data$sample_id, tsne$sample_id),]
  tsne_desc = data.frame(base = "all",
                         name = "All Cells UMAP")
  
  write_feather(tsne, file.path(mappingFolder, "tsne.feather"))
  write_feather(as_tibble(tsne_desc), file.path(mappingFolder, "tsne_desc.feather"))
}
