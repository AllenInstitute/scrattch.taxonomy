#' Read in a reference data set in Allen taxonomy format
#'
#' @param taxonomyDir Directory containing the AIT file -OR- a direct h5ad file name -OR- a URL of a publicly accessible AIT file.
#' @param anndata_file If taxonomyDir is a directory, anndata_file must be the file name of the anndata object (.h5ad) to be loaded in that directory. If taxonomyDir is a file name or a URL, then anndata_file is ignored.
#' @param log.file.path Path to write log file to. Defaults to current working directory. 
#' @param read_subsampled Should a subsampled version of the taxonomy be loaded (default=FALSE). If TRUE, additional parameters for read_subsampled_h5ad must be provided.
#' @param ... Additional parameters passed to read_subsampled_h5ad (note that "h5ad" is inferred from anndata_file above)
#'
#' @return Organized reference object ready for mapping against.
#' 
#' @import anndata
#'
#' @export
loadTaxonomy = function(taxonomyDir = getwd(), 
                        anndata_file = "AI_taxonomy.h5ad",
                        log.file.path=getwd(),
                        force=FALSE,
                        read_subsampled=FALSE,
                        ...){
  
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
      options(timeout = 9000) 
      download.file(anndata_file, file.path(taxonomyDir, anndata_object), timeout=9000)
    }
    anndata_file = anndata_object
  }
  
  ## Load from directory name input 
  if(file.exists(file.path(taxonomyDir, anndata_file))){
    if(read_subsampled){
      print("Loading SUBSAMPLED reference taxonomy into memory from .h5ad")
      ## Load subsampled taxonomy directly!
      h5ad = file.path(taxonomyDir, anndata_file)
      AIT.anndata = read_subsampled_h5ad(h5ad,...)
    } else {
      print("Loading COMPLETE reference taxonomy into memory from .h5ad")
      ## Load taxonomy directly!
      AIT.anndata = read_h5ad(file.path(taxonomyDir, anndata_file))
    }
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
  AIT.anndata$uns$title <- gsub(".h5ad","",anndata_file)

  ## Return
  return(AIT.anndata)
}





#' Read an h5ad file into R after subsampling in python
#'
#' @param h5ad The location of an h5ad file (if not local, it will be downloaded).
#' @param cluster_vector Optional subsampling spec. Can be one of:
#'   (1) character vector of cell IDs to keep,
#'   (2) logical vector TRUE=keep/FALSE=omit (aligned to obs order),
#'   (3) numeric vector of 1-based row indices to keep.
#'   If provided, cluster_column and subsample are ignored.
#' @param cluster_column The column name in obs where clusters are stored
#'   (default = "cluster_id"). Ignored if cluster_vector is provided.
#' @param subsample How much subsampling should be done before running statistics.
#'   Default is 100 per item in cluster_column. Ignored if cluster_vector is provided.
#' @param seed The seed to use for reproducibility.
#' @param python_exe Python executable (default = "python3").
#' @param keep_files Should the subsampled h5ad be kept on disk?
#' @param outdir Output directory for temporary or kept files.
#' @param read_backed If TRUE, read subsampled h5ad with backed="r".
#'
#' @return the subsampled AnnData object
#'
#' @import anndata
#' @export
read_subsampled_h5ad <- function(
    h5ad,
    cluster_vector = NULL,
    cluster_column = "cluster_id",
    subsample = 100,
    seed = 42,
    python_exe = "python3",
    keep_files = FALSE,
    outdir = ifelse(keep_files, getwd(), tempdir()),
    read_backed = FALSE
) {
  stopifnot(is.character(h5ad), length(h5ad) == 1)
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  
  # ------------------------------------------------------------
  # 0) Ensure local h5ad
  # ------------------------------------------------------------
  if (grepl("^https?://", h5ad)) {
    local_h5ad <- file.path(outdir, basename(h5ad))
    if (!file.exists(local_h5ad)) {
      options(timeout = 9000)
      download.file(h5ad, destfile = local_h5ad, mode = "wb", timeout = 9000)
    }
  } else {
    local_h5ad <- h5ad
  }
  if (!file.exists(local_h5ad))
    stop("h5ad file not found: ", local_h5ad)
  
  # ------------------------------------------------------------
  # 1) Determine subsample IDs in R
  # ------------------------------------------------------------
  keep_ids <- NULL
  
  if (is.null(cluster_vector)) {
    ## ---------- Mode 1: compute subsample internally ----------
    if (!requireNamespace("anndata", quietly = TRUE))
      stop("R package 'anndata' is required.")
    
    if (!exists("subsampleCells", mode = "function"))
      stop("subsampleCells() not found in R environment.")
    
    ad <- anndata::read_h5ad(local_h5ad, backed = "r")
    obs <- ad$obs
    if (!(cluster_column %in% colnames(obs)))
      stop("cluster_column not found in obs: ", cluster_column)
    
    cl <- obs[, cluster_column]
    names(cl) <- rownames(obs)
    
    keep_mask <- subsampleCells(cl, subSamp = subsample, seed = seed)
    if (!is.logical(keep_mask) || length(keep_mask) != nrow(obs))
      stop("subsampleCells() must return a logical vector aligned to obs.")
    
    keep_ids <- names(cl)[keep_mask]
    
  } else {
    ## ---------- Mode 2: subsample provided in advance ----------
    if (is.character(cluster_vector)) {
      keep_ids <- as.character(cluster_vector)
      
    } else if (is.logical(cluster_vector)) {
      ad <- anndata::read_h5ad(local_h5ad, backed = "r")
      if (length(cluster_vector) != nrow(ad$obs))
        stop("Logical cluster_vector length must equal number of cells.")
      keep_ids <- rownames(ad$obs)[cluster_vector]
      
    } else if (is.numeric(cluster_vector)) {
      ad <- anndata::read_h5ad(local_h5ad, backed = "r")
      keep_ids <- rownames(ad$obs)[as.integer(cluster_vector)]
      
    } else {
      stop("cluster_vector must be character, logical, or numeric.")
    }
  }
  
  if (length(keep_ids) == 0)
    stop("No cells selected for subsampling.")
  
  # ------------------------------------------------------------
  # 2) Write subsample_ids.csv
  # ------------------------------------------------------------
  ids_csv <- file.path(
    outdir,
    sprintf("subsample_ids_%d.csv", length(keep_ids))
  )
  write.csv(
    data.frame(cell_id = keep_ids),
    ids_csv,
    row.names = FALSE,
    quote = TRUE
  )
  
  # ------------------------------------------------------------
  # 3) Python helper: subset h5ad by IDs
  #    - Reads non-backed (avoids X-access errors in backed mode)
  #    - Forces sp.csc_matrix (not csc_array) + float64 for X & raw
  #      so R reads dgCMatrix with double-typed x slot
  # ------------------------------------------------------------
  py_script <- file.path(outdir, "subsample_h5ad.py")
  writeLines(c(
    "#!/usr/bin/env python3",
    "import argparse, numpy as np, pandas as pd",
    "import scipy.sparse as sp, anndata as ad",
    "",
    "def csc64(X):",
    "  '''Guarantee scipy.sparse.csc_matrix with float64 dtype.'''",
    "  if X is None: return None",
    "  m = sp.csc_matrix(X) if sp.issparse(X) else sp.csc_matrix(np.asarray(X))",
    "  m = m.astype(np.float64)",
    "  m.sort_indices()",
    "  return m",
    "",
    "def main():",
    "  ap = argparse.ArgumentParser()",
    "  ap.add_argument('--h5ad', required=True)",
    "  ap.add_argument('--ids',  required=True)",
    "  ap.add_argument('--out',  required=True)",
    "  a = ap.parse_args()",
    "  ids = pd.read_csv(a.ids)['cell_id'].astype(str).tolist()",
    "  seen = set(); ids = [x for x in ids if not (x in seen or seen.add(x))]",
    "  adata = ad.read_h5ad(a.h5ad)          # full read; avoids backed-mode X issues",
    "  sub = adata[ids, :].copy()",
    "  sub.X = csc64(sub.X)",
    "  if sub.raw is not None:",
    "    r = sub.raw.to_adata()",
    "    raw = ad.AnnData(X=csc64(r.X), obs=sub.obs.copy(), var=r.var.copy())",
    "    raw.obs_names = sub.obs_names.copy()",
    "    raw.var_names = r.var_names.copy()",
    "    sub.raw = raw",
    "  for k in list(sub.layers.keys()):",
    "    sub.layers[k] = csc64(sub.layers[k])",
    "  sub.write_h5ad(a.out, compression='gzip')",
    "  print(f'[done] {len(ids)} cells  X={type(sub.X).__name__}  dtype={sub.X.dtype}')",
    "",
    "if __name__ == '__main__': main()"
  ), py_script)
  
  out_h5ad <- file.path(outdir, paste0("subsampled_", basename(local_h5ad)))
  if (!grepl("\\.h5ad$", out_h5ad, ignore.case = TRUE))
    out_h5ad <- paste0(out_h5ad, ".h5ad")
  
  res <- system2(
    python_exe,
    args = c(
      py_script,
      "--h5ad", local_h5ad,
      "--ids", ids_csv,
      "--out", out_h5ad
    ),
    stdout = TRUE,
    stderr = TRUE
  )
  if (!file.exists(out_h5ad))
    stop("Python subsampling failed:\n", paste(res, collapse = "\n"))
  
  # ------------------------------------------------------------
  # 4) Read subsampled h5ad back into R
  # ------------------------------------------------------------
  ad_sub <- if (isTRUE(read_backed)) {
    anndata::read_h5ad(out_h5ad, backed = "r")
  } else {
    anndata::read_h5ad(out_h5ad)
  }
  
  if (!keep_files) {
    suppressWarnings(file.remove(py_script))
    suppressWarnings(file.remove(ids_csv))
  }
  
  ad_sub
}

