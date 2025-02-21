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

## Provide hierarchy of the taxonomy
## -- This will be used for all mapping algorithms unless otherwise specified
## -- This MUST be from broadest to most specific types, and NOT vice versa
hierarchy = list("broad_type" = 1, "primary_type_label" = 2)

## Compute top 1000 binary marker genes for clusters (or use a pre-existing vector)
binary.genes = top_binary_genes(taxonomy.counts, taxonomy.anno$primary_type_label, 1000)

## Compute UMAP coordinates (or use precomputed coordinates)
pcs  <- prcomp(logCPM(taxonomy.counts)[binary.genes,], scale = TRUE)$rotation
umap.coords = umap(pcs[,1:30])$layout

## Set rownames to your annotation and UMAP data.frames with sample identifiers (Required!)
rownames(taxonomy.anno) = taxonomy.anno$sample_name
rownames(umap.coords) = colnames(taxonomy.counts)

## Identify Ensembl IDs 
# Common NCBI taxIDs: Human = 9606; Mouse = 10090; Macaque (rhesus) = 9544; Marmoset = 9483
ensembl_id <- geneSymbolToEnsembl(gene.symbols = rownames(taxonomy.counts), ncbi.taxid = 10090)

## Align taxonomy metadata with AIT standard
# NOTE: for this particular example, nothing gets updated
full.taxonomy.anno <- updateTaxonomyMetadata(taxonomy.anno)
taxonomy.anno      <- full.taxonomy.anno$metadata

## Build Allen Insitute Taxonomy, for large taxonomies you can pass in tpm and cluster_stats if pre-computed.
AIT.anndata = buildTaxonomy(meta.data = taxonomy.anno,
                            title = "Tasic2016",
                            counts = as(taxonomy.counts, "dgCMatrix"),
                            normalized.expr = NULL,
                            highly_variable_genes = NULL,
                            marker_genes = list("marker_genes_binary" = binary.genes),
                            ensembl_id = ensembl_id,
                            cluster_stats = NULL, ## Pre-computed cluster stats
                            embeddings = list("X_umap" = umap.coords),
                            ##
                            dend = NULL, ## Pre-computed dendrogram
                            taxonomyDir = getwd(), ## This is where our taxonomy will be created
                            hierarchy = hierarchy,
                            ##
                            subsample=2000)

## Check whether the taxonomy file is valid
AIT.anndata$uns$valid = checkTaxonomy(AIT.anndata)


## Create Shiny directory (AIBS-internal)
createShiny(AIT.anndata,
            shinyDir = getwd(),  # Replace location with desired location for shiny directory output
            metadata_names = NULL)
```
