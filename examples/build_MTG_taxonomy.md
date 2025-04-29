# Tutorial: Building an AIT taxonomy of human MTG 

In this tutorial we demonstrate how to setup an Allen Institute Taxonomy object using scrattch.taxonomy--see the [schema definitions](https://github.com/AllenInstitute/AllenInstituteTaxonomy/tree/main/schema) for more details. In this tutorial we will use data from middle temporal gyrus of give neurotypical human donors, along with associated metadata and cell type assignments created as part of the Sealle Alzheimer's disease brain cell atlas (SEA-AD). These data are publicly accessible at https://portal.brain-map.org/atlases-and-data/rnaseq/human-mtg-10x_sea-ad. We then create a child taxonomy only including neuronal cells for patch-seq mapping (in scrattch.mapping and scrattch.taxonomy R libraries).

These data are already QCed and nicely packaged in h5ad (counts and metadata) and an associated dendrogram files. "cluster_label", "subclass_label", and "class_label" correspond to SEA-AD supertype, subclass, and class, respectively, and are used for defining the hierarchy.  

*We strongly encourage running this code within the scrattch docker environment.  This example was created using docker://alleninst/scrattch:1.1.2 and will likely fail if run using any earlier scrattch versions.*

#### Prepare taxonomy data set:

First we download the data, read it into R, and subset the data set to 100 cells per cluster (to reduce computational burden and more evenly sample cell types).

```R
## Load scrattch.taxonomy
library(scrattch.taxonomy)
library(reticulate)
cell_type_mapper <- import("cell_type_mapper")

## Define and go to your working directory
taxonomyDir = getwd() # Replace with location of taxonomy
if(!file.exists(taxonomyDir)) dir.create(taxonomyDir)
setwd(taxonomyDir)

## Download the reference data to the working directory and read it in
seaad_url  <- "https://sea-ad-single-cell-profiling.s3.us-west-2.amazonaws.com/MTG/RNAseq/Reference_MTG_RNAseq_final-nuclei.2022-06-07.h5ad"
dend_url   <- "https://brainmapportal-live-4cc80a57cd6e400d854-f7fdcae.divio-media.net/filer_public/0f/37/0f3755cb-3acb-4b93-8a62-5d6adc74c673/dend.rds"
#download.file(seaad_url,"Reference_MTG_RNAseq_final-nuclei.2022-06-07.h5ad")  # NOTE: we recommend downloading via the web browser, as this command may fail
#download.file(dend_url,"Reference_MTG_dend.rds")
seaad_data <- read_h5ad("Reference_MTG_RNAseq_final-nuclei.2022-06-07.h5ad")
seaad_dend <- readRDS("Reference_MTG_dend.rds")

## Subsample data (this can be done either here, or within buildTaxonomy)
keepCells <- subsampleCells(seaad_data$obs$cluster_label,100,seed=42)

## Get (subsampled) subset data and annotations
taxonomy.counts = (seaad_data$X)[keepCells,]
cn <- c("sample_name","cluster_label","cluster_confidence","subclass_label","class_label",
        "external_donor_name_label","age_label","donor_sex_label")
taxonomy.metadata = seaad_data$obs[keepCells,cn]

## Ensure count matrix and annotations are in the same order (this shouldn't be needed)
taxonomy.metadata = taxonomy.metadata[match(rownames(taxonomy.counts), taxonomy.metadata$sample_name),]
colnames(taxonomy.metadata) <- gsub("_label","",colnames(taxonomy.metadata))
```

#### Align to AIT schema

Next, we update the metadata fields to align with the AIT schema.

```R
## Set up the levels of hierarchy for all mapping functions later
## -- This MUST be from broadest to most specific types, and NOT vice versa
## -- Note that "cluster" is the SEAAD supertypes and will be named "cluster_id" below
hierarchy = list("class", "subclass", "cluster_id")

## Identify Ensembl IDs 
# Common NCBI taxIDs: Human = 9606; Mouse = 10090; Macaque (rhesus) = 9544; Marmoset = 9483
ensembl_id <- geneSymbolToEnsembl(gene.symbols = colnames(taxonomy.counts), ncbi.taxid = 9606)

## Update the metadata to align with AIT schema
colnames(taxonomy.metadata)[colnames(taxonomy.metadata)=="cluster"]             = "cluster_id"
colnames(taxonomy.metadata)[colnames(taxonomy.metadata)=="donor_sex"]           = "self_reported_sex"
colnames(taxonomy.metadata)[colnames(taxonomy.metadata)=="external_donor_name"] = "donor_id"
taxonomy.metadata$load_id           = "Not reported"
taxonomy.metadata$assay             = "10x 3' v3"  
taxonomy.metadata$organism          = "Homo sapiens"
taxonomy.metadata$anatomical_region = "Middle temporal gyrus"
taxonomy.metadata$suspension_type   = "nucleus"
taxonomy.metadata$is_primary_data   = TRUE
taxonomy.metadata$self_reported_ethnicity = "unknown"

## Check taxonomy metadata aligngs with AIT standard and perform minor error corrections
## Also add ontology terms corresponding to the above schema elements (and can also correct misspellings, etc.)
full.taxonomy.anno <- computeOntologyTerms(taxonomy.metadata, standardize.metadata=TRUE, print.messages=TRUE) 
# NOTE: We encourage reviewing the messages from this function CAREFULLY, as some assumptions are made when calculating ontology terms

## Save final metadata data frame
taxonomy.anno <- full.taxonomy.anno$metadata
```


#### Build parent taxonomy:

By default the taxonomy will now be set up for mapping algorithms.

```R
## Build Allen Insitute Taxonomy, for large taxonomies you can pass in tpm and cluster_stats if pre-computed.
AIT.anndata = buildTaxonomy(title="SEAAD_MTG",
                            meta.data = taxonomy.anno,
                            hierarchy = hierarchy,
                            ## --- Optional parameters ---
                            counts = taxonomy.counts,
                            normalized.expr = NULL,
                            highly_variable_genes = 1000, ## Select top 1000 binary genes
                            marker_genes = NULL,
                            ensembl_id = ensembl_id,
                            cluster_stats = NULL, ## Pre-computed cluster stats
                            embeddings = "highly_variable_genes_standard", # Compute UMAP coordinates internally
                            ##
                            dend = seaad_dend, ## Pre-computed dendrogram
                            taxonomyDir = getwd(), ## This is where our taxonomy will be created
                            addMapMyCells = TRUE, 
                            ##
                            add.dendrogram.markers = TRUE,  # Allow tree mapping. Very slow, but required for downstream patch-seq analysis.
                            subsample=100)

## Check whether the taxonomy file is valid (This happens within buildTaxonomy and is not strictly necessary)
AIT.anndata = checkTaxonomy(AIT.anndata)
```


#### Build child taxonomy (with neurons):

Now let's create a version of the taxonomy which is compatible with patchseqQC and can be filtered to remove off target cells from mapping. In this case we are going to exclude all the non-neuronal cells since we know that (1) patch-seq cells are all neurons and (2) the transcriptomics of patch-seq cells are often contaminated by non-neuronal transcripts that can "trick" algorithms into thinking the cells are non-neuronal. 

```R
## Identify the neuronal clusters
neuron.cells = AIT.anndata$obs$class!="Non-neuronal and Non-neural"

## Build the mode surrounding these neurons
AIT.anndata = buildTaxonomyMode(AIT.anndata, 
                                mode.name = "neurons", 
                                highly_variable_genes = 1000,
                                embeddings = "highly_variable_genes_neurons",
                                retain.cells = neuron.cells, 
                                subsample = 100, 
                                overwrite = TRUE)
```

