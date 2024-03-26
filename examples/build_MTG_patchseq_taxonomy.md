
# Tutorial: Building a Patch-seq Shiny taxonomy - Human MTG

In this example we demonstrate how to setup a Patch-seq Shiny taxonomy using scrattch.mapping for viewing on MolGen Shiny and running mapping algorithms against. This tutorial parallels the other tutorial for building a Patch-seq Shiny taxonomy, but using the [Hodge et al taxonomy](https://www.nature.com/articles/s41586-019-1506-7) as reference and query Patch-seq data from [Berg et al 2020](https://www.nature.com/articles/s41586-021-03813-8). If you are bringing your own data to the tutorial you can replace all sections that say "**FOR EXAMPLE ONLY**" with your own data munging steps. 

For creating a standard taxonomy and mapping against it, the following input variables are required: 

To create a taxonomy in AIT format:
* Reference count matrix (gene x cell), with genes as rownames and sample identifiers as colnames
* Reference annotation data.frame (cell x field), with sample identifiers as rownames (this should include cluster calls)
* Variable and/or marker genes (vector), can be calculated if not provided
* UMAP coordinates (cell x 2), can be calculated if not provided
* Dendrogram of clusters (dendrogram), optional (except for tree mapping), can be calculated if not provided

To map patch-seq (or other) data to the taxonomy:
* Query count matrix and metadata (of the format above)
* Query annotation data.frame (cell x field), optional, will be passed through to MolGen Shiny directory

## Part 1: Building the taxonomy

This section describes how to build a reference taxonomy for mapping. If you have already built a taxonomy in AIT format, you can skip to Part 2: Mapping to the taxonomy.

### 1.1: Read in the reference data

The first step is to read in your cell (columns) by gene (rows) data matrix of counts alongside the metadata associated with each cell.  These are saved as taxonomy.counts and taxonomy.metadata, respectively in the section below.  We use data from Hodge et al 2019 **FOR EXAMPLE ONLY**.  Note that colnames(taxonomy.counts) and rownames(taxonomy.metadata) should be identical and correspond to cell IDs, and rownames(taxonomy.counts) should be the gene names (in the is case, gene symbols).

```R
## Load the library
suppressPackageStartupMessages({
  library(data.table) # for using fread below
  library(hodge2019data) # devtools::install_github("AllenInstitute/hodge2019data")
  library(R.utils) # install.packages('R.utils')'  # for using fread below
}) 

## Load the complete dendrogram from this paper (saved in hodge2010data)
data(dend_Hodge2019) 

## Download data and metadata from the website (this is slow)
## NOT AVAILABLE YET. Instead copy Hodge2019_metadata.csv and Hodge2019_counts.csv.gz from "/allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/Taxonomies/AIT15.3/" to your working directory.

## Read data and metadata into R
taxonomy.metadata <- read.csv("/allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/Taxonomies/AIT30/Hodge2019_metadata.csv",row.names=1)
taxonomy.counts   <- as.matrix(fread("/allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/Taxonomies/AIT30/Hodge2019_counts.csv.gz"),rownames=1) ## This requires R.utils
colnames(taxonomy.counts) <- rownames(taxonomy.metadata) # To correct "-" to "." conversion introduced at some point. 
```

### 1.2: Create the (parent) AIT Taxonomy 
=======
### Create the base Shiny Taxonomy for the ENTIRE Hodge et al 2019 data set

```R
## This is where our taxonomy will be created
taxonomy = "/allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/Taxonomies/AIT15.3/"
```

This section will create the parent taxonomy for the reference data.  In this case, we include up 1000 cells for **every** cell type defined in Hodge et al 2019, along with their associated metadata, and will subsample the clusters and cells further at a later step.

This code block loads the scrattch.taxonomy library, and then calculates variables genes and defines a UMAP **using a very basic approach**.  If variable genes and/or 2-dimensional coordinates already exist, they can be provided to buildTaxonomy below rather than calculated in this way. 

```R
#devtools::install_github("AllenInstitute/scrattch.taxonomy", ref = "KL_div")
suppressPackageStartupMessages({
  library(scrattch.taxonomy)
})

## Compute top 1000 binary marker genes for clusters
binary.genes = top_binary_genes(taxonomy.counts, taxonomy.metadata$cluster_label, 1000)

## Compute UMAP coordinates
pcs  <- prcomp(logCPM(taxonomy.counts)[binary.genes,], scale = TRUE)$rotation
umap.coords = umap(pcs[,1:30])$layout

## Set rownames to your annotation and UMAP data.frames (Required!)
rownames(umap.coords) = colnames(taxonomy.counts)
```

The next step builds the parent taxonomy using a single call to the function buildTaxonomy.  After running this script, your taxonomy will contain all of the data and metadata in standard formats, and will be ready for correlation and tree mapping.  However, **you still need to run buildPatchseqMapping in the next section** for tree mapping, subsetting, and other QC metrics.  

```R
## This is where our taxonomy will be created
taxonomy = "/allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/Taxonomies/AIT30/"

## Build Shiny taxonomy 
AIT.anndata = buildTaxonomy(counts = taxonomy.counts,
              meta.data   = taxonomy.metadata,
              dend        = dend_Hodge2019,  # If this is omitted buildTaxonomy will generate a dendrogram
              feature.set = binary.genes,
              umap.coords = umap.coords,
              taxonomyName= "MTG_Hodge2019",
              taxonomyDir = taxonomy,
              subsample   = 1000) # Minimal subsampling (N=1000)

## Alternatively, if you have already created the taxonomy, you can load it using "loadTaxonomy"
## Load the taxonomy (from h5ad file name)
#AIT.anndata = loadTaxonomy(taxonomy, "MTG_Hodge2019.h5ad")
```

### 1.3: Create a child taxonomy for patch-seq mapping

Now let's create a version of the taxonomy which is compatible with patchseqQC and can be filtered to remove off target cells from mapping. **You are creating a new version of the base taxonomy which can be reused by specifying the provided `mode.name` in `scrattch.taxonomy::mappingMode()` as discussed next.**

```R
## Setup the taxonomy for patchseqQC to infer off.target contamination
AIT.anndata = buildPatchseqTaxonomy(AIT.anndata,
                                    mode.name = "AIT30.1", ## Give a name to off.target filtered taxonomy
                                    subsample = 100, ## Subsampling for the new taxonomy.
                                    subclass.column = "subclass_label", 
                                    class.column = "class_label", ## The column by which off-target types are determined.
                                    off.target.types = c("Non-neuronal"), ## The off-target class.column labels for patchseqQC.
                                    subclass.subsample = 100, ## Subsampling is for PatchseqQC contamination calculation.
                                    num.markers = 50, ## Number of markers for each annotation in `class_label`
                                    taxonomyDir = taxonomy) # This will create a subfolder in the reference taxonomy directory
```

The `buildPatchseqTaxonomy` function does the following, updating the anndata variable and file accordingly:
* Creates a "child" taxonomy that only includes the cells and clusters requested (e.g., neurons for patch-seq mapping)
* Defines marker genes for each node of the dendogram for use with tree mapping
* Creates a subsetted version of the parent dendrogram that also includes node marker genes
* Creates the required marker and expression variables for 'QC_markers' for use with patchseqQC
* Creates a table of cell to cluster probability (e.g., 'membership') values for calculation of KL divergence and creation of constellation diagrams
Currently a subfolder for the child taxonomy is also created, but everything in that folder is also stored in the anndata, and so it can be safely ignored.

**At this point the reference taxonomy is created and ready for Patch-seq mapping.**

## Part 2: Mapping to the taxonomy

The rest of this example demonstrates how to map patch-seq data to the reference defined above.  This is a primary use case of the scrattch.mapping library, and is entirely distinct from creation of a reference taxonomy in Part 1 above. 

### 2.1: Read in the QUERY (e.g., Patch-seq) data

The first step is to read in your cell (columns) by gene (rows) data matrix of **log-normalized** query counts alongside (optional) metadata associated with each cell.  These are saved as query.logCPM and query.anno, respectively in the section below.  We use data from Berg et al 2020 **FOR EXAMPLE ONLY**.  Note that colnames(query.logCPM) and rownames(query.anno) should be identical and correspond to cell IDs, and rownames(query.logCPM) should be the gene names (in the is case, gene symbols).  **Only common gene names in query.logCPM and taxonomy.counts will be used for mapping.**


```R
## Load scrattch.mapping
suppressPackageStartupMessages({
  library(scrattch.mapping)
})

## Download data and metadata from GitHub repo for Berg et al 2022
download.file("https://github.com/AllenInstitute/patchseq_human_L23/raw/master/data/input_patchseq_data_sets.RData", "patchseq.RData", mode="wb")

## Load the data
load("patchseq.RData")

## Rename the query data and metadata for convenience
query.anno = annoPatch  # Some cell annotations for all cells from the paper 
query.logCPM = datPatch # logCPM values for all cells from the paper
```


### 2.2: Set scrattch.mapping mode

**Do not skip this step!** Here we set taxonomy mode to the relevant "child" taxonomy defined in 1.3 above.  This will let the mapping functions know to  only consider the unfiltered cells and cell types (e.g., use subsetted cells from neuronal clusters for Patch-seq mapping).

```R
AIT.anndata = mappingMode(AIT.anndata, mode="AIT30.1")
```

### 2.3: Map query cells to Patch-seq reference

This is the step that performs the mapping.  Since we have very different query and reference data sets, we are omitting Seurat mapping in this example. Note that each mapping algorithm can map to multiple levels of the taxonomy. 

```R
# This function is part of the 'scrattch.mapping' library
query.mapping = taxonomy_mapping(AIT.anndata= AIT.anndata,
                                 query.data = query.logCPM, 
                                 corr.map   = TRUE, # Flags for which mapping algorithms to run
                                 tree.map   = TRUE, 
                                 seurat.map = FALSE, 
                                 label.cols = c("cluster_label", "subclass_label" ,"class_label")) # Columns to map against from AIT.anndata$obs
```

### 2.4 Prepare Patch-seq Shiny directory

This section performs additional query data QC and then data outputs the files necessary for visualization of Patch-seq data with molgen-shiny tools into a folder. More specifically, the function buildMappingDirectory:
* Creates a new folder to deposit all relevant query files (this is what you copy into the the molgen-shiny directory)
* Performs patchSeqQC on the query data to define Normalized Marker Sum (NMS) and other QC scores (see https://github.com/PavlidisLab/patchSeqQC for details)
* Calculates KL Divergence and associated Core/I1/I2/I3/PoorQ calls used to help assess Patch-seq quality
* Calculates and outputs tree mapping probabilities of every cell mapping to every cluster and tree node
* Calculates and outputs UMAP coordinates for query cells integrated into reference UMAP space
*Note for users from outside the Allen Institute: a set of molgen-shiny visualization tools are only available at the Allen Institute.  However, all of the above variables are output as .feather files that can be read into R using the function read_feather. Please submit issues if more information is needed.*

```R
buildMappingDirectory(AIT.anndata    = AIT.anndata, 
                      mappingFolder  = "/allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/Taxonomies/AIT30/TEST",
                      query.data     = query.logCPM, ## Counts or CPM are required here, but function can convert from log to linear values
                      query.metadata = query.anno,
                      query.mapping  = query.mapping,
                      doPatchseqQC   = TRUE)  ## Set to FALSE if not needed or if buildPatchseqTaxonomy was not run.
```
