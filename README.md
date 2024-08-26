# scrattch.taxonomy

Generalized taxonomy building scripts for RNA-seq based taxonomies following the [Allen Institute schema](https://github.com/AllenInstitute/scrattch.taxonomy/tree/main/schema).

**A list of available taxonomies in this format is available at [Taxonomy_list.md](https://github.com/AllenInstitute/scrattch.taxonomy/blob/main/Taxonomy_list.md)**. As of 10 November 2023, these represent published single-cell RNAseq taxonomies from the Allen Institute and other groups, and are largely complementary to taxonomies included at the [Brain Knowledge Platform](https://portal.brain-map.org/atlases-and-data/bkp), although efforts to integrate are ongoing.

## Documentation

You can find a detail description of all scrattch.taxonomy functions here: ![Documentation](https://github.com/AllenInstitute/scrattch.taxonomy/blob/main/scrattch.taxonomy_0.1.pdf)

Update notes are here: ![Versions](https://github.com/AllenInstitute/scrattch.taxonomy/blob/main/VERSIONS.md)

## Docker

We have setup a docker environemnt for scrattch.taxonomy and scrattch.mapping that contains all the required dependencies and the current version of both scrattch.taxonomy and scrattch.mapping. This docker is accessible through docker hub via: `njjai/scrattch_mapping:0.6.2`.

#### HPC usage:

##### Non-interactive
`singularity shell --cleanenv docker://njjai/scrattch_mapping:0.6.3 Rscript YOUR_CODE.R`

##### Interactive
`singularity shell --cleanenv docker://njjai/scrattch_mapping:0.6.3`


## Installation

While we advise using the provided docker, you can also install scrattch.taxonomy directly from GitHub as follows:

```
# Quickly, but without the vignettes:
devtools::install_github("AllenInstitute/scrattch.taxonomy")

# More slowly, but with the vignettes:
devtools::install_github("AllenInstitute/scrattch.taxonomy", build_vignettes=TRUE, force=TRUE)
```

This strategy **might not work outside the docker** due to complicated dependencies. Also note that `doMC` may need to be installed manually from the download link at https://r-forge.r-project.org/R/?group_id=947 if you use Windows. Vignettes are provided below.

## Usage examples

1. [**Build a Shiny taxonomy**](https://github.com/AllenInstitute/scrattch.taxonomy/blob/main/examples/build_taxonomy.md) This example provides the basics for creating a new taxonomy compatible with scrattch.mapping mapping functions and (internal Allen Institute) MolGen Shiny tools.
   
## Reporting issues

If you run into any issues, please let Nelson and Jeremy know or [**create a new issue in the 'Issues' tab above**](https://github.com/AllenInstitute/scrattch-taxonomy/issues).

## TODO

- [ ] Update documentation.

## Done
