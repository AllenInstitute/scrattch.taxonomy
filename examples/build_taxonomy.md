# Tutorial: Building a Shiny taxonomy 

In this tutorial we demonstrate how to setup a Shiny taxonomy using scrattch.taxonomy for running mapping algorithms against with scrattch.mapping and for viewing on (internal Allen Institute) MolGen Shiny tools. 

#### Required inputs:

* Count matrix (gene x cell), with genes as rownames and sample identifiers as colnames
* Annotation data.frame (cell x field), with sample identifiers as rownames
* Variable and/or marker genes (vector), can be calculated if not provided
* UMAP coordinates (cell x 2), can be calculated if not provided

#### Additional prerequisites:

* Installation of the `tasic2016data` data package [from here](https://github.com/AllenInstitute/tasic2016data/), or replace with your own data set.

#### Build taxonomy:

```R
## Load scrattch.taxonomy
library(scrattch.taxonomy)

## Load in example count data and annotations (or replace with your own)
library(tasic2016data)
taxonomy.counts = tasic_2016_counts
taxonomy.anno = tasic_2016_anno

## Ensure count matrix and annotations are in the same order.
taxonomy.anno = taxonomy.anno[match(colnames(taxonomy.counts), taxonomy.anno$sample_name),]

## Ensure 'cluster' field exists, as required by scrattch.taxonomy.
taxonomy.anno$cluster = taxonomy.anno$broad_type

## Compute top 1000 binary marker genes for clusters (or use a pre-existing vector)
binary.genes = top_binary_genes(taxonomy.counts, taxonomy.anno$cluster, 1000)

## Compute UMAP coordinates (or use precomputed coordinates)
pcs  <- prcomp(logCPM(taxonomy.counts)[binary.genes,], scale = TRUE)$rotation
umap.coords = umap(pcs[,1:30])$layout

## Set rownames to your annotation and UMAP data.frames with sample identifiers (Required!)
rownames(taxonomy.anno) = taxonomy.anno$sample_name
rownames(umap.coords) = colnames(taxonomy.counts)

## This is where our taxonomy will be created
# NOTE: replace 'taxonomyDir' location below with desired output folder location
taxonomyDir = "/allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/10x_seq/tasic_2016/"

## Build Shiny taxonomy 
AIT.anndata = buildTaxonomy(counts = taxonomy.counts,
                meta.data = taxonomy.anno,
                feature.set = binary.genes,
                umap.coords = umap.coords,
                taxonomyName = "Tasic2016", ## NEW!
                taxonomyDir = taxonomyDir,
                subsample=2000)

## Add markers to dendrogram
AIT.anndata = addDendrogramMarkers(AIT.anndata = AIT.anndata)

```

# Setup MapMyCells taxonomy

In this tutorial we demonstrate how to setup a MapMyCells taxonomy using scrattch.taxonomy for running HANN mapping algorithm against with scrattch.mapping.

#### Required inputs:

* Hierarchy of the taxonomy as a list, such as class_label, subclass_label, cluster_label
* AIT.anndata, MolGen Shiny (scrattch) taxonomy

#### Additional prerequisites:

* Installation of the `tasic2016data` data package [from here](https://github.com/AllenInstitute/tasic2016data/), or replace with your own data set.

#### Build taxonomy:

```R
library(reticulate)
cell_type_mapper <- import("cell_type_mapper")

## Provide hierarchy of the taxonomy
hierarchy = list("broad_type_label", "primary_type_label")

## Build MapMyCells taxonomy
AIT.anndata = buildHANNMapMyCells(AIT.anndata, hierarchy)
```
