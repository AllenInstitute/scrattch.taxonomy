# scrattch.taxonomy

Generalized taxonomy building scripts for RNA-seq based taxonomies following the [Allen Institute schema](https://github.com/AllenInstitute/AllenInstituteTaxonomy).

## Documentation

You can find a detail description of all scrattch.taxonomy functions here: ![Documentation](https://github.com/AllenInstitute/scrattch-taxonomy/blob/main/scrattch.taxonomy_0.1.pdf)

Update notes are here: ![Versions](https://github.com/AllenInstitute/scrattch-taxonomy/blob/main/VERSIONS.md)

## Docker

We have setup a docker environemnt for scattch.taxonomy and scrattch.mapping that contains all the required dependencies and the current version of both scrattch.taxonomy and scrattch.mapping. This docker is accessible through docker hub via: `bicore/scrattch_mapping:latest`.

#### HPC usage:

##### Non-interactive
`singularity exec --cleanenv docker://bicore/scrattch_mapping:latest Rscript YOUR_CODE.R`

##### Interactive
`singularity shell --cleanenv docker://bicore/scrattch_mapping:latest`


## Installation

While we advice using the provided docker, you can also install scrattch.taxonomy directly from github as follows:

*Note: slight edits to installation will be needed while repo is private.  Also note that `doMC` may need to be installed manually from the download at https://r-forge.r-project.org/R/?group_id=947 if you use Windows.*

```
# Quickly, but without the vignettes:
devtools::install_github("AllenInstitute/scrattch-taxonomy")

# More slowly, but with the vignettes:
devtools::install_github("AllenInstitute/scrattch-taxonomy", build_vignettes=TRUE, force=TRUE)
```

Note that this strategy might not work outside the docker due to complicated dependencies. Vignettes are provided below.

## Usage examples

1. [**Build a Shiny taxonomy**](https://github.com/AllenInstitute/scrattch-taxonomy/blob/main/examples/build_taxonomy.md) This examples provides the basics for creating a new Shiny taxonomy compatible with MolGen shiny and scrattch.mapping.

2. [**Build and map against a small mouse PatchSeq taxonomy**](https://github.com/AllenInstitute/scrattch-taxonomy/blob/main/examples/build_patchseq_taxonomy.md) This examples provides the basics for updating a Shiny taxonomy to be compatible with patchseq style mapping and visualization on MolGen Shiny.

3. [**Build and map against a human MTG PatchSeq taxonomy**](https://github.com/AllenInstitute/scrattch-taxonomy/blob/main/examples/build_MTG_patchseq_taxonomy.md) This examples provides a simplified example for creating a standard taxonomy and updating it to be compatible with patchseq style mapping and visualization on MolGen Shiny.  Data is from Hodge et al for human MTG and patch-seq examples are from Berg et al (layer 2-3, excitatory neurons). 

## Library vignettes

*Note: links below will not work while repo is private.*

1. [**Build reference directory for mapping.**](http://htmlpreview.github.io/?https://github.com/AllenInstitute/mfishtools/blob/master/vignettes/build_reference_taxonomy.html)  This vignette provides an example of how to convert a *completed* single cell RNA-seq analysis (e.g., a counts matrix + cell type assignments) into a standard reference taxonomy. Resulting taxonomy files are used as input for various mapping techniques in this package, and are also compatible with tools for visualiation of taxonomies at the Allen Institute and Patch-seq QC and visualization. **This process must be run first.**  

2. [**Map patch-seq data and output directory.**](http://htmlpreview.github.io/?https://github.com/AllenInstitute/mfishtools/blob/master/vignettes/complete_patchseq_analysis.html)  This vignette goes through how to map a small data set against a reference taxonomy. Here we use a subset of tasic2016data as an example but the intention is for mapping of patch-seq data.  

## Reporting issues

If you run into any issues, please let Nelson and Jeremy know or [**create a new issue in the 'Issues' tab above**](https://github.com/AllenInstitute/scrattch-taxonomy/issues).

## TODO

- [ ] Update documentation.

## Done
