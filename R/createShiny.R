#' This function builds the minimum files required for Shiny
#'
#' Precomputed clusters must be provided.  In the anndata object these will be stored using the term "cluster".  If celltypeColumn is anything other than cluster, then any existing "cluster" column will be overwritten by celltypeColumn.  Values can be provided without colors and ids (e.g., "cluster") or with them (e.g., "cluster_label" + "cluster_color" + "cluster_id").  In this case cluster_colors is ignored and colors are taken directly from the metadata.  Cluster_id's will be overwritten to match dendrogram order.
#'
#' @param AIT.anndata A reference taxonomy anndata object.
#' @param shinyDir The location to save files for MolGen Shiny (AIBS)
#' @param metadata_names An optional named character vector where the vector NAMES correspond to columns in the metadata matrix and the vector VALUES correspond to how these metadata should be displayed in Shiny. This is used for writing the desc.feather file later.
#' 
#' @import scrattch.hicat 
#' @import feather
#' @import tibble
#' @import dplyr
#' @import Matrix
#' @import pvclust
#' @import anndata
#'
#' @return AIT anndata object in the specified format (only if return.anndata=TRUE)
#'
#' @export
createShiny = function(AIT.anndata,
                        shinyDir,
                        metadata_names = NULL){

  ## Create shinyDir if needed
  if(!file.exists(shinyDir))
    dir.create(shinyDir)
  
  ## Pull filtering / subsampling for current mode
  kpSub = !AIT.anndata$uns$filter[[AIT.anndata$uns$mode]]

  ## Get the data and metadata matrices
  counts = as(Matrix::t(AIT.anndata[kpSub,]$raw[["X"]]),"dgCMatrix")
  colnames(counts) = rownames(AIT.anndata$obs)
  rownames(counts) = rownames(AIT.anndata$var)
  meta.data = AIT.anndata[kpSub,]$obs
  umap.coords = AIT.anndata[kpSub,]$obsm[[AIT.anndata$uns$default_embedding]];   # THIS NEEDS TO BE UPDATED when we update default_embedding to a list
  rownames(umap.coords) = rownames(meta.data)

  ## Create log TPM matrix from counts
  if(!is.null(AIT.anndata$X)){
    tpm.matrix = Matrix::t(AIT.anndata[kpSub,]$X)
  }else{
    tpm.matrix = logCPM(counts)
  }

  ## If needed, set "cluster_id" as the first level of the hierarchy
  hier = as.data.frame(AIT.anndata$uns$hierarchy)
  meta.data$cluster_id = meta.data[,colnames(hier)[which(hier==max(hier))]]
  
  ## Next, generate and output the data matrices. Make sure matrix class is dgCMatrix and not dgTMatrix.  
  sample_id <- meta.data$sample_id <- colnames(tpm.matrix)
  gene      = rownames(tpm.matrix)
  cluster   = meta.data$cluster_id   
  names(cluster) = rownames(meta.data) ## For get_cl_medians

  ## Counts feather - Not used for shiny, but useful for saving raw data nonetheless
  print("===== Building counts.feather =====")
  counts.tibble = as_tibble(counts)
  counts.tibble = cbind(gene, counts.tibble)
  counts.tibble = as_tibble(counts.tibble)
  feather::write_feather(counts.tibble, file.path(shinyDir, "counts.feather"))

  ## Data_t feather
  print("===== Building data_t.feather =====")
  data.tibble = as_tibble(tpm.matrix)
  data.tibble = cbind(gene, data.tibble)
  data.tibble = as_tibble(data.tibble)
  feather::write_feather(data.tibble, file.path(shinyDir, "data_t.feather"))
  
  ## Data feather
  print("===== Building data.feather =====")
  norm.data.t = t(as.matrix(tpm.matrix))
  norm.data.t = as_tibble(norm.data.t)
  norm.data.t = cbind(sample_id, norm.data.t)
  norm.data.t = as_tibble(norm.data.t)
  feather::write_feather(norm.data.t, file.path(shinyDir, "data.feather"))
  
  ## Output tree
  dend = json_to_dend(AIT.anndata$uns$dend[[AIT.anndata$uns$mode]])  # Might need to deal with modes here.  For now we write the dendrogram for the current mode
  saveRDS(dend, file.path(shinyDir, "dend.RData"))

  ## Output tree order
  outDend = data.frame(cluster = labels(dend), order = 1:length(labels(dend)))
  write.csv(outDend, file.path(shinyDir, "ordered_clusters.csv"))
    
  ## Minor reformatting of metadata file, then write metadata file
  print("===== Format the metadata =====")
  meta.data$sample_id = colnames(counts)
  colnames(meta.data)[colnames(meta.data)=="sample_name"] <- "sample_name_old" # Rename "sample_name" to avoid shiny crashing
  colnames(meta.data)[colnames(meta.data)=="cluster_id"]  <- "cluster" 
  meta.data = .formatMetadata(meta.data, cluster_colors, "sample_id")
  meta.data$cluster_id <- as.numeric(factor(meta.data$cluster_label,levels=labels(dend))) # Reorder cluster ids to match dendrogram
  feather::write_feather(meta.data, file.path(shinyDir,"anno.feather"))

  ## Write the desc file.  
  anno_desc = create_desc(meta.data, use_label_columns = TRUE)
  # Subset the desc file to match metadata_names, if provided
  if(!is.null(metadata_names)){
    desc <- anno_desc[match(names(metadata_names), as.character(as.matrix(anno_desc[,1]))),]
    desc[,2] <- as.character(metadata_names[as.character(as.matrix(desc[,1]))])
    desc <- desc[!is.na(desc$base),]  # Remove missing values
    anno_desc <- desc
  }
  feather::write_feather(anno_desc, file.path(shinyDir,"desc.feather"))
  
  ## Write the UMAP coordinates.  
  print("===== Building umap/tsne feathers (precalculated) =====")
  tsne      = data.frame(sample_id = rownames(umap.coords),
                         all_x = umap.coords[,1],
                         all_y = umap.coords[,2])
  tsne      = tsne[match(meta.data$sample_id, tsne$sample_id),]
  tsne_desc = data.frame(base = "all",
                         name = paste("Embedding:",AIT.anndata$uns$default_embedding)) # THIS NEEDS TO BE UPDATED when we update default_embedding to a list

  write_feather(tsne, file.path(shinyDir,"tsne.feather"))
  write_feather(tsne_desc, file.path(shinyDir,"tsne_desc.feather"))

  ##
  print("===== Building count, median, sum feathers =====")
  all_clusters = unique(meta.data$cluster_label) ## cluster_id
  count_gt0 = matrix(0, ncol = length(all_clusters), nrow = nrow(data.tibble))
  count_gt1 = sums = medianmat = count_gt0

  ## Compute the number of genes with greater than 0,1 counts, gene sums and medians
  for (i in 1:length(all_clusters)) {
    cluster = all_clusters[i]
    cluster_samples = which(meta.data$cluster_label == cluster)
    cluster_data    = tpm.matrix[,cluster_samples,drop=F]
    cluster_counts  = counts[,cluster_samples,drop=F]
    count_gt0[, i]  = Matrix::rowSums(cluster_counts > 0)
    count_gt1[, i]  = Matrix::rowSums(cluster_counts > 1)
    sums[, i]       = Matrix::rowSums(cluster_counts)
    medianmat[, i]  = apply(cluster_data, 1, median)
  }

  ##
  all_cl_label <- meta.data$cluster_id[match(all_clusters,meta.data$cluster_label)]
  all_cl_label <- paste0("cluster_",all_cl_label) # to get format in 'cluster_[cluster_id]'
  colnames(count_gt0) = colnames(count_gt1) = colnames(sums) = colnames(medianmat) = all_cl_label # all_clusters ## allClust
  count_gt0 = cbind(gene = gene, as.data.frame(count_gt0))
  count_gt1 = cbind(gene = gene, as.data.frame(count_gt1))
  sums = cbind(gene = gene, as.data.frame(sums))
  medianmat = cbind(gene = gene, as.data.frame(medianmat))

  ##
  count_n = meta.data %>% 
              arrange(cluster_id) %>% 
              group_by(cluster_id) %>% 
              summarise(n_cells = n())

  ##
  feather::write_feather(count_gt0, file.path(shinyDir,"count_gt0.feather"))
  feather::write_feather(count_gt1, file.path(shinyDir,"count_gt1.feather"))
  feather::write_feather(count_n,   file.path(shinyDir,"count_n.feather"))
  feather::write_feather(medianmat, file.path(shinyDir,"medians.feather"))
  feather::write_feather(sums,      file.path(shinyDir,"sums.feather"))
}



  # ## ----------
  # ## Run auto_annotate, this changes cell_id to cell_id.
  # meta.data$cell_id = colnames(counts)
  # meta.data = .formatMetadata(meta.data, cluster_colors)

  #   ## Gather clusters
  # clustersUse = unique(meta.data$cluster_label)
