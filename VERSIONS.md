## scrattch.taxonomy v1.1.2

Update to allow counts to be provided as cellxgene or genexcell.

### Major changes
* Creation of logCPM_byRows to allow counts to be input as cellxgene or genexcell and to avoid multiple large matrix transpositions
* Update loadTaxonomy to auto-calulate normalized counts so these don't need to be saved in h5ad files
* Addition of a gene.meta.data slot in buildTaxonomy for gene information

### Minor changes
* A few additional bug fixes from v1.1.1
* Minor updates to the MTG example 
* Output additional warnings when metadata is changed with updateTaxonomyMetadata


--

## scrattch.taxonomy v1.1.1

Update to better integrate the Allen Institute schema and downstream scrattch.mapping and scrattch.patchseq functionality.

### Major changes
* Creation of buildTaxonomyMode function for adding flexible modes/subsets/children (without adding patchseq QC)
* Update dendrogram and MapMyCells handling with schema
* Allow computing of binary genes and UMAP within buildTaxonomy and buildTaxonomyMode
* Creation of a new example using buildTaxonomyMode with human MTG

### Minor changes
* Bug fixes
* Speed up for geneSymbolToEnsembl function
* Update of brain region atlases to correct ontology links
* Update self_reported_ethnicity_ontology_term_id to deal with "unknown" and "multiethnic" entries
* Minor changes to error handling 

--

## scrattch.taxonomy v0.9.1

Major update to integrate the new Allen Institute schema

### Major changes
* Introduction of AIT schema into code
* Many new functions to deal with schema, ontologies, taxonomy checks, etc.
* Updated example to use new schema
* Metadata error correction


### Minor changes
Bug fixes

--

## scrattch.taxonomy v0.7.1

Updates to improve correlation mapping

### Major changes
* New function `updateHighlyVariableGenes` to update genes used in mapping

### Minor changes
* Synchronized versions between scrattch.taxonomy, scrattch.mapping, an scrattch.patchseq
* Set mode-specific variable genes for "standard"
* Minor (but breaking) bug fix for AIT.anndata <---> AIT_anndata
* Updated documentation

--

## scrattch.taxonomy v0.5.14

Upgrades to streamline process

### Major changes
* Updated example
* Pull code for setting up tree and hierarchical mapping into buildTaxonomy

### Minor changes
Bug fixes

--

## scrattch.taxonomy v0.5.13

Supporting changes for scrattch.mapping updates to hierarchical and Seurat mapping

### Major changes
* Update docker to downgrade Seurat and SeuratObject packages to older versions
* Building the hierarchical mapping statistics so it can be correctly stored for multiple modes  

### Minor changes
Bug fix allowing users to provide their own dendrogram

--

## scrattch.taxonomy v0.2

## Major changes

* Inclusion of the taxonomy schema on the GitHub page here.
* Removal of vignettes and integration into website examples
* Inclusion of additional functions for usability on Windows
* Moved dendrogram from Hodge et al to new data package (https://github.com/AllenInstitute/hodge2019data)

## Minor updates

* Update examples to reflect splitting of scrattch.taxonomy and scrattch.mapping libraries
* Small changes for typos and error correction
