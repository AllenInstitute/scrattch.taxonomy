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
