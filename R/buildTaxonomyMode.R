#' Save marker genes for patchSeqQC
#'
#' This function saves all the variables required for applying the patchseq QC algorithm `pathseqtools` (which is an more flexible version of the `patchSeqQC` algorithm) to AIT.anndata$uns. This is only used for patch-seq analysis.  Requirements for input include:
# ----- Subclass calls for each cell
# ----- Broad class class calls for each cell
# ----- Distinction of neuron (e.g., mappable type) vs. non-neuron (e.g., contamination type)
#'
#' @param AIT.anndata A reference taxonomy anndata object.
#' @param subsample The number of cells to retain per cluster (default = 100).
#' @param subclass.column Column name corresponding to the moderate-resolution cell types used for the cell types of interest (default = "subclass_label").
#' @param class.column Column name corresponding to the low-resolution cell types used for the off-target cell types (default = "class_label").
#' @param off.target.types A character vector of off-target (also known as 'contamination') cell types.  This must include at least one of the cell types found in "class.column" and/or "subclass.column" (both columns are checked)
#' @param mode.name A name to identify the new taxonomy version.
#' @param subclass.subsample The number of cells to retain for PatchseqQC contamination calculation (default = 100, probably no need to change).
#' @param num.markers The maximum number of markers to calculate per node per direction (default = 50)
#' @param taxonomyDir The location to save shiny output (default = current working directory).
#' @param ... Additional variables to be passed to `addDendrogramMarkers`
#' 
#' The following variables are added to AIT.anndata$uns:  
#' $dend[[mode.name]]  
#' $filter[[mode.name]]  
#' $QC_markers[[mode.name]]  
#' ...$markers,  
#' ...$countsQC,  
#' ...$cpmQC,  
#' ...$classBr,  
#' ...$subclassF,  
#' ...$allMarkers  
#' ...$de_genes  
#' $memb[[mode.name]]  
#' ...$memb.ref,  
#' ...$map.df.ref  
#' 
#' @import patchseqtools
#' @import scrattch.hicat
#'
#' @return AIT.anndata An updated AIT.anndata variable with the above content added to AIT.anndata$uns for the relevant mode.name.
#'
#' @export
buildTaxonomyMode = function(AIT.anndata,
                                 ...){

  ## Allow user to create modes which reduce the number of cells or genes.
  ## This is useful for creating a taxonomy that can be used for patchseq or spatial mapping.
  ## We are calling these children taxonomy of the full parent.

}
