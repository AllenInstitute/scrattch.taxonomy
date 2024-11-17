# Tutorial: Building a Shiny taxonomy 

In this tutorial we demonstrate how to setup a Shiny taxonomy using scrattch.taxonomy for running mapping algorithms against with scrattch.mapping and for viewing on (internal Allen Institute) MolGen Shiny tools. It uses prepackaged data from [Tasic et al 2016](https://www.nature.com/articles/nn.4216)

#### Required inputs:

* Count matrix (gene x cell), with genes as rownames and sample identifiers as colnames
* Annotation data.frame (cell x field), with sample identifiers as rownames
* Variable and/or marker genes (vector), can be calculated if not provided
* UMAP coordinates (cell x 2), can be calculated if not provided

#### Additional prerequisites:

* Installation of the `tasic2016data` data package [from here](https://github.com/AllenInstitute/tasic2016data/), or replace with your own data set.

#### Build taxonomy:

By default the taxonomy will now be set up for all four mapping algorithms

```R
## Load scrattch.taxonomy
library(scrattch.taxonomy)
library(reticulate)
cell_type_mapper <- import("cell_type_mapper")

## Load in example count data and annotations (or replace with your own)
library(tasic2016data)
taxonomy.counts = tasic_2016_counts
taxonomy.anno   = tasic_2016_anno
taxonomy.anno   = taxonomy.anno[match(colnames(taxonomy.counts), taxonomy.anno$sample_name),]
keep            = taxonomy.anno$broad_type!="Unclassified"
taxonomy.counts = taxonomy.counts[,keep]
taxonomy.anno   = taxonomy.anno[keep,]

## Ensure 'cluster' field exists, as required by scrattch.taxonomy.
## -- Note that this is called "cluster" here, but will be called "cluster_label" everywhere else
taxonomy.anno$cluster = taxonomy.anno$primary_type_label

## Provide hierarchy of the taxonomy
## -- This will be used for all mapping algorithms unless otherwise specified
## -- This MUST be from broadest to most specific types, and NOT vice versa
hierarchy = list("broad_type_label", "cluster_label")

## Compute top 1000 binary marker genes for clusters (or use a pre-existing vector)
binary.genes = top_binary_genes(taxonomy.counts, taxonomy.anno$cluster, 1000)

## Compute UMAP coordinates (or use precomputed coordinates)
pcs  <- prcomp(logCPM(taxonomy.counts)[binary.genes,], scale = TRUE)$rotation
umap.coords = umap(pcs[,1:30])$layout

## Set rownames to your annotation and UMAP data.frames with sample identifiers (Required!)
rownames(taxonomy.anno) = taxonomy.anno$sample_name
rownames(umap.coords) = colnames(taxonomy.counts)

## This is where our taxonomy will be created
## -- NOTE: replace 'taxonomyDir' location below with desired output folder location
taxonomyDir = getwd()

## Build Allen Insitute Taxonomy, for large taxonomies you can pass in tpm and cluster_stats if pre-computed.
AIT.anndata = buildTaxonomy(counts = as(taxonomy.counts, "dgCMatrix"),
                            tpm = NULL,
                            meta.data = taxonomy.anno,
                            feature.set = binary.genes,
                            umap.coords = umap.coords,
                            taxonomyDir = taxonomyDir,
                            taxonomyTitle = "Tasic2016",
                            hierarchy = hierarchy,
                            subsample=200)

## Create Shiny directory (AIBS-internal)
createShiny(AIT.anndata,
            shinyDir = getwd(),  # Replace location with desired location for shiny directory output
            metadata_names = NULL)
```
