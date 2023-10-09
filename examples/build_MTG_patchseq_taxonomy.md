
# Tutorial: Building a Patch-seq Shiny taxonomy - Human MTG

In this example we demonstrate how to setup a Patch-seq Shiny taxonomy using scrattch.mapping for viewing on MolGen Shiny and running mapping algorithms against. This tutorial parallels the other tutorial for building a Patch-seq Shiny taxonomy, but using the [Hodge et al taxonomy](https://www.nature.com/articles/s41586-019-1506-7) as reference and query Patch-seq data from [Berg et al 2020](https://www.nature.com/articles/s41586-021-03813-8).  

## Part 1: Building the taxonomy

#### Required inputs:

* Standard Shiny taxonomy setup following the "build_taxonomy" tutorial (repeated below)
* Query Patch-seq count matrix and metadata (example provided below)

### Read in the REFERENCE MTG information (Hodge et al 2019)

```R
## Load the library
library(scrattch.taxonomy)
library(data.table) # for using fread below

## Load the complete dendrogram from this paper (saved in scrattch.taxonomy)
data(dend_Hodge2019) 

## Download data and metadata from the website (this is slow)
## NOT AVAILABLE YET. Instead copy Hodge2019_metadata.csv and Hodge2019_counts.csv.gz from "/allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/Taxonomies/AIT15.3/" to your working directory.

## Read data and metadata into R
taxonomy.metadata <- read.csv("/allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/Taxonomies/AIT15.3/Hodge2019_metadata.csv",row.names=1)
taxonomy.counts   <- as.matrix(fread("/allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/Taxonomies/AIT15.3/Hodge2019_counts.csv.gz"),rownames=1) ## This requires R.utils
colnames(taxonomy.counts) <- rownames(taxonomy.metadata) # To correct "-" to "." conversion introduced at some point. 
```

### Create the base Shiny Taxonomy for the ENTIRE Hodge et al 2019 data set

```R
## This is where our taxonomy will be created
taxonomy = "/allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/Taxonomies/AIT15.3/"

For creating the standard taxonomy (see "Building a Shiny taxonomy" tutorial):
* Count matrix (gene x cell), with genes as rownames and sample identifiers as colnames
* Annotation data.frame (cell x field), with sample identifiers as rownames
* Variable and/or marker genes (vector), can be calculated if not provided
* UMAP coordinates (cell x 2), can be calculated if not provided

For mapping of patch-seq data:
* Query patchseq count matrix and metadata (of the format above; example provided below)

#### Additional prerequisites:

* Installation of the `hodge2019data` data package [from here](https://github.com/AllenInstitute/hodge2019data/), or replace with your own data set.
* Installation of `scrattch.mapping` [from here](https://github.com/AllenInstitute/scrattch.mapping) for data mapping. 

#### Read in the REFERENCE MTG information (Hodge et al 2019)
```R
## Load the libraries
library(scrattch.taxonomy)
library(scrattch.mapping)
library(hodge2019data)

## Read data and metadata into R
taxonomy.metadata <- metadata_Hodge2019
taxonomy.counts   <- data_Hodge2019
# dend_Hodge2019 is also part of this hodge2019data and will be used below
```

#### Create the base Shiny Taxonomy for the full Hodge et al 2019 data set (neurons + non-neurons)
```R
## This is where our taxonomy will be created
# Replace this location below with the folder location you'd like to sue
taxonomy = "MTG_Hodge2019/"

## Compute top 1000 binary marker genes for clusters
binary.genes = top_binary_genes(taxonomy.counts, taxonomy.metadata$cluster_label, 1000)

## Compute UMAP coordinates
pcs  <- prcomp(logCPM(taxonomy.counts)[binary.genes,], scale = TRUE)$rotation
umap.coords = umap(pcs[,1:30])$layout

## Set rownames to your annotation and UMAP data.frames (Required!)
rownames(umap.coords) = colnames(taxonomy.counts)

## Build Shiny taxonomy 
AIT.anndata = buildTaxonomy(counts = taxonomy.counts,
              meta.data   = taxonomy.metadata,
              dend        = dend_Hodge2019,  # If this is omitted buildTaxonomy will generate a dendrogram
              feature.set = binary.genes,
              umap.coords = umap.coords,
              taxonomyName= "MTG_Hodge2019",
              taxonomyDir = taxonomy,
              subsample   = 250
## Load the taxonomy (from h5ad file name)
AIT.anndata = loadTaxonomy(taxonomy, "AI_taxonomy.h5ad")
```

### Build the patchseq taxonomy

Now let's create a version of the taxonomy which is compatible with patchseqQC and can be filtered to remove off target cells from mapping. **You are creating a new version of the base taxonomy which can be reused by specifying the provided `mode.name` in `scrattch.taxonomy::mappingMode()` as dicusssed next.**

```R
## Setup the taxonomy for patchseqQC to infer off.target contamination
AIT.anndata = buildPatchseqTaxonomy(AIT.anndata,
                                    mode.name = "patchseq", ## Give a name to off.target filterd taxonomy
                                    subsample = 100, ## Subsampling is only for PatchseqQC contamination calculation.
                                    subclass.column = "subclass_label", 
                                    class.column = "class_label", ## The column by which off-target types are determined.
                                    off.target.types = c("Non-neuronal"), ## The off-target class.column labels for patchseqQC.
                                    num.markers = 50, ## Number of markers for each annotation in `class_label`
                                    taxonomyDir = taxonomy) # This will create a subfolder in the reference taxonomy directory
```
The `buildPatchseqTaxonomy` function return/created the following:

* An updated AIT.anndata object for Patch-seq mapping and QC steps.
* Created the required marker and expression variables for 'QC_markers' and save under 'mode.name' in the uns
* Created the required cell to cluster 'membership' variables and save under 'mode.name' in the uns
* An updated dendrogram and a saved in the 'mode.name' subdirectory

Now let's check to make sure the anndata object is formatted correctly.  This can manually be done using the function "checkTaxonomy" (as below) but also happens automatically above.

```R
checkTaxonomy(AIT.anndata,taxonomyDir)
```
You can check the log file if directed, but if the value returned is "TRUE", then the taxonomy should work for downstream scrattch.taxonomy and scrattch.mapping functions and **at this point the reference taxonomy is created and ready for Patch-seq mapping.**

## Part 2: Mapping to the taxonomy

The rest of this example demonstrates how to read in Patch-seq data (using [Berg et al 2020](https://www.nature.com/articles/s41586-021-03813-8)) and map it to the neuronal cell types from the [Hodge et al 2019 taxonomy](https://www.nature.com/articles/s41586-019-1506-7).  Note that we need to load the scrattch.mapping library for this part. 

### First read in and process QUERY Patch-seq data

```R
## Load scrattch.mapping
library(scrattch.mapping)

## Download data and metadata from GitHub repo for Berg et al 2022
download.file("https://github.com/AllenInstitute/patchseq_human_L23/raw/master/data/input_patchseq_data_sets.RData", "patchseq.RData", mode="wb")

## Load the data
load("patchseq.RData")

## Rename the query data and metadata for convenience
query.anno = annoPatch  # Some cell annotations for all cells from the paper 
query.logCPM = datPatch # logCPM values for all cells from the paper
```


#### Set scrattch.mapping mode

Now we will set scrattch.mapping to use only cells not in the off.target.types, this will filter the taxonomy and adjust the dendrogram to remove any `cluster` in `off.target.types`.
```R
AIT.anndata = mappingMode(AIT.anndata, mode="patchseq")
```

#### Map against the patchseq taxonomy

This is the step that performs the mapping.  Since we have very different query and reference data sets, we are omitting Seurat mapping in this example. 

```R
# This function is part of the 'scrattch.mapping' library
query.mapping = taxonomy_mapping(AIT.anndata= AIT.anndata,
                                 query.data = query.logCPM, 
                                 corr.map   = TRUE, # Flags for which mapping algorithms to run
                                 tree.map   = TRUE, 
                                 seurat.map = FALSE, 
                                 label.cols = c("cluster_label", "subclass_label" ,"class_label")) # Columns to map against from AIT.anndata$obs
```

### Setup the Patch-seq Shiny taxonomy files for human MTG:

This step outputs the files necessary for visualization of Patch-seq data with molgen-shiny tools.  *Note that this code block also generates additional QC metrics including NMS with PatchseqQC.  To perform patch-seq QC without building this directory, use the `applyPatchseqQC` function.*

```R
buildMappingDirectory(AIT.anndata    = AIT.anndata, 
                      mappingFolder  = "/allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/Taxonomies/AIT15.3/TEST",
                      query.data     = query.logCPM, ## Counts or CPM are required here, but function can convert from log to linear values
                      query.metadata = query.anno,
                      query.mapping  = query.mapping,
                      doPatchseqQC   = TRUE)  ## Set to FALSE if not needed or if buildPatchseqTaxonomy was not run.
```
