#' Add marker genes to reference dendrogram for tree mapping
#'
#' @param AIT.anndata A reference taxonomy anndata object.
#' @param mode.name Taxonomy mode to determine which version of filtering to use.
#' @param celltypeColumn Column name correspond to the cell type names in the dendrogram. At least two cells per cell type in the dendrogram must be included. Default is to take the highest resolution level of the hierarchy
#' @param subsample The number of cells to retain per cluster (default = 100)
#' @param num.markers The maximum number of markers to calculate per pairwise differential calculation per direction (default = 20)
#' @param de.param Differential expression (DE) parameters for genes and clusters used to define marker genes.  By default the values are set to the 10x nuclei defaults from scrattch.hicat, except with min.cells=2 (see notes below).
#' @param calculate.de.genes Default=TRUE. If set to false, the function will search for "de_genes" in a file called "de_genes_[mode].RData" in the folder taxonomyDir, and will use those if available.
#' @param save.shiny.output Should standard output files be generated and saved to the AIT file (default=TRUE).  These are not required for tree mapping, but are required for building a patch-seq shiny instance. See notes.
#' @param mc.cores Number of cores to use for running this function to speed things up.  Default = 1.  Values>1 are only supported in an UNIX environment and require `foreach` and `doParallel` R libraries.
#' @param bs.num Number of bootstrap runs for creating the dendrogram (default of 100)
#' @param p proportion of marker genes to include in each iteration of the mapping algorithm.
#' @param low.th the minimum difference in Pearson correlation required to decide on which branch
#' @param overwriteMarkers If markers already are calculated a tree, should they be overwritten (default = TRUE)
#' @param taxonomyDir The location to create the directory with taxonomy information. Both the taxonomy and the de.genes_[mode.name].RData get saved here.
#'
#' By default VERY loose parameters are set for de_param in an effort to get extra marker genes for each node.  The defaults previously proposed for 10x nuclei are the following `de_param(low.th = 1, padj.th = 0.01, lfc.th = 1, q1.th = 0.3, q2.th = NULL, q.diff.th = 0.7, de.score.th = 100, min.cells = 2, min.genes = 5)`. See the function `de_param` in the scrattch.hicat for more details.  
#'
#' If save.shiny.output=TRUE, the following two values will get saved to the uns:
#'       `memb.ref`   - matrix indicating how much confusion there is the mapping between each cell all of the nodes in the tree (including all cell types) when comparing clustering and mapping results with various subsamplings of the data
#'       `map.df.ref` - Result of tree mapping for each cell in the reference against the clustering tree, including various statistics and marker gene evidence.  This is the same output that comes from tree mapping.
#' 
#' @import feather
#' @import scrattch.hicat
#' @import MatrixGenerics
#' @import dendextend
#'
#' @return An updated dendrogram variable that is the same as `dend` except with marker genes added to each node.
#'
#' @export
addDendrogramMarkers = function(AIT.anndata,
                                mode.name = NULL,
                                mode = AIT.anndata$uns$mode, # For back compatibility. Will be used if mode.name is not provided
                                celltypeColumn = NULL, 
                                subsample = 100,
                                num.markers = 20,
                                de.param=scrattch.hicat::de_param(low.th = 0.1,
                                                                  padj.th = 1,
                                                                  lfc.th = 0.01,
                                                                  q1.th = 0.1,
                                                                  q2.th = NULL,
                                                                  q.diff.th = 0.1,
                                                                  de.score.th = 1,
                                                                  min.cells = 1,
                                                                  min.genes = 1),
                                calculate.de.genes = TRUE,
                                save.shiny.output = TRUE,
                                mc.cores=1, 
                                bs.num=100, 
                                p=0.8, 
                                low.th=0.1,
                                overwriteMarkers = TRUE,
                                taxonomyDir = file.path(AIT.anndata$uns$taxonomyDir),
                                ...){  # Not used, but required to include
  suppressWarnings({ # wrapping the whole function in suppressWarnings to avoid having this printed a zillion times: 'useNames = NA is deprecated. Instead, specify either useNames = TRUE or useNames = FALSE.'
  
  ## Deal with back compatibility of mode
  if(is.null(mode.name)) mode.name = mode
  rm(mode)

  ## File name for de_genes_file
  de_genes_file <- file.path(taxonomyDir,paste0("de_genes_",mode.name,".RData"))

  ## Create the taxonomyDir if somehow it doesn't exist yet
  if(!dir.exists(taxonomyDir)){ dir.create(taxonomyDir) }

  ## Determine the cluster column if not provided
  if(is.null(celltypeColumn)){
    hierarchy = AIT.anndata$uns$hierarchy
    if(is.null(hierarchy)) {stop("Either a cell type column or a hierarchy must be provided to addDendrogramMarkers")}
    hierarchy = hierarchy[order(unlist(hierarchy))]
    if(is.null(hierarchy)) stop("Hierarchy must be included in the standard AIT mode in proper format to create a mode.  Please run checkTaxonomy().")
    celltypeColumn = names(hierarchy)[length(hierarchy)][[1]]
  }
  
  ## If counts are included but normalized counts are not, calculate normalized counts
  if((!is.null(AIT.anndata$raw$X))&(is.null(AIT.anndata$X))){
    normalized.expr = log2CPM_byRow(AIT.anndata$raw$X)
    AIT.anndata$X   = normalized.expr 
  }
  
  ##
  if(is.null(AIT.anndata$X)){ stop("No data found in AIT.anndata$X. The full count matrix is required for determining marker genes.") }

  ## Filter and Subsample
  keep.samples = ((!AIT.anndata$uns$filter[[mode.name]]) & subsampleCells(AIT.anndata$obs[[celltypeColumn]], subsample)) ##  & is.element(cluster.vector, labels(dend))

  ## Checks and data formatting
  dend = json_to_dend(AIT.anndata$uns$dend[[mode.name]])

  ## norm.data
  norm.data = Matrix::t(AIT.anndata$X[AIT.anndata$obs_names[keep.samples],])
  
  ## metadata 
  metadata = AIT.anndata$obs[keep.samples,] %>% as.data.frame()

  ## celltype labels
  cluster.vector = setNames(metadata[,celltypeColumn], rownames(metadata))
  cluster.vector = factor(cluster.vector, levels=labels(dend))
  
  ## Define and collect marker genes
  print("Define some relevant variables")
  cl.df = as.data.frame(metadata[match(labels(dend), cluster.vector),])
  cl.df$cluster_label = cl.df[,celltypeColumn]
  rownames(cl.df) = 1:length(labels(dend)) 
  cl.label  = as.factor(setNames(cl.df$cluster_label, rownames(cl.df)))
  select.cl = droplevels(as.factor(setNames(match(metadata[,celltypeColumn], cl.label), AIT.anndata$obs_names[keep.samples])))
  
  ## CHECK IF THIS IS NEEDED - CONFIRMED. THIS WILL BREAK IF NOT INCLUDED
  ## We might need to relabel the dendrogram from 1 to #clusters in order
  labels(dend) = names(cl.label)[match(labels(dend), cl.label)]
  
  ## Compute markers
  print("Define marker genes and gene scores for the tree")
  if((sum(!is.na(get_nodes_attr(dend, "markers"))) == 0) | (overwriteMarkers == TRUE)){
    if(calculate.de.genes|(!file.exists(de_genes_file))){
      if(calculate.de.genes&(!file.exists(de_genes_file)))
        print(paste(de_genes_file,"does not exist.  Need to define marker genes!"))
      print("=== NOTE: This step can be very slow (several minute to many hours).")
      print("      To speed up the calculation (or if it crashes) try decreasing the value of subsample.")
      de.genes = scrattch.hicat::select_markers(norm.dat=norm.data, 
                                cl=select.cl, 
                                n.markers= num.markers, 
                                de.param = de.param, 
                                de.genes = NULL)$de.genes

      save(de.genes, file=de_genes_file)
    } else {
      ## Check that file exists before loading
      print("...using precomputed marker genes.")
      load(de_genes_file)
    }

    ## Check number of markers for each leaf
    min.marker.gene.count = as.numeric(as.character(lapply(de.genes, function(x) x$num)))
    if((sum(min.marker.gene.count<2,na.rm=TRUE)>0)|(sum(is.na(min.marker.gene.count)>0)))({
      stop("Marker genes could not be calculated for at least one node in the tree. Tree mapping will not work in this situation. We recommend loosening the de.param parameters and trying again, but warn that some cell types may not be well resolved.")
    })

    ## Note: this is not a robust way to address this issue!
    print("Calculate gene scores")
    gene.score = scrattch.hicat::get_gene_score(de.genes)
  
    print("Build the reference dendrogram")
    invisible(capture.output({  # Avoid printing lots of numbers to the screen
      reference = build_reference(cl=select.cl, 
                                  norm.dat=norm.data, 
                                  dend=dend, 
                                  de.genes=de.genes, 
                                  cl.label=cl.label, 
                                  up.gene.score=gene.score$up.gene.score, 
                                  down.gene.score=gene.score$down.gene.score, 
                                  n.markers=num.markers)
    }))
    labels(reference$dend) = setNames(colnames(reference$cl.dat), colnames(reference$cl.dat))
    print("...marker gene calculation for reference complete")
  } else {
    ## Use already defined marker genes on tree
    print("Build the reference dendrogram")
    reference <- list(
      dend = dend,
      cl.dat = get_cl_means(norm.data, select.cl)
    )
    print("...existing markers used in reference")
  }
    
  ##
  if(sum(!is.na(get_nodes_attr(reference$dend, "original_label")) > 0)){
    print("This section is needed if the starting dendrogram is from the nomenclature GitHub ")
    reference$dend = revert_dend_label(reference$dend,get_nodes_attr(reference$dend, "original_label"),"label")
  }

  #print("Save the reference dendrogram for this mode")  # NOTE: we don't need to save this dendrogram, as it is saved within the AIT file
  #dend = reference$dend
  #saveRDS(dend, file.path(taxonomyModeDir, "dend.RData"))
  ## Note, this overwrites the initial dendrogram but has slightly different formatting from the read, which could cause issues
  ## dend = readRDS(AIT.anndata$uns$dend[[mode.name]])
    
  ##
  if(save.shiny.output){
    ## NOTE: These are used for KL mapping and potentially for future constellation diagrams
    print("Build membership table of reference vs. reference for use with patch-seq mapping")
    invisible(capture.output({  # Avoid printing lots of numbers to the screen
      memb.ref   = .map_dend_membership(reference$dend, 
                                        reference$cl.dat, 
                                        map.dat=as.matrix(norm.data), 
                                        map.cells=names(select.cl),
                                        mc.cores=mc.cores, 
                                        bs.num=bs.num, 
                                        p=p, 
                                        low.th=low.th)
      map.df.ref = summarize_cl(reference$dend, 
                                memb.ref, 
                                as.matrix(norm.data))
    }))
    memb.ref   = memb.ref[AIT.anndata$obs_names[keep.samples],]
    map.df.ref = map.df.ref[AIT.anndata$obs_names[keep.samples],]
    
    AIT.anndata$uns$memb[[mode.name]]$memb.ref = as.data.frame.matrix(memb.ref)
    AIT.anndata$uns$memb[[mode.name]]$map.df.ref = map.df.ref
    #save(memb.ref, map.df.ref, file=file.path(taxonomyModeDir, "membership_information_reference.rda")) # No need to save this to a file either
  }

  ##
  print("Save the dendrogram into .h5ad")
  AIT.anndata$uns$dend[[mode.name]] = toJSON(dend_to_json(reference$dend))
  AIT.anndata$uns$title <- gsub(".h5ad","",AIT.anndata$uns$title)
  AIT.anndata$write_h5ad(file.path(AIT.anndata$uns$taxonomyDir, paste0(AIT.anndata$uns$title, ".h5ad")))
  
  ##
  return(AIT.anndata)
  
}) # End suppressWarnings

}
