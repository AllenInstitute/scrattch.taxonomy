# scrattch.taxonomy

Generalized taxonomy building scripts for RNA-seq based taxonomies following the [Allen Institute schema](https://github.com/AllenInstitute/scrattch.taxonomy/tree/main/schema).

**A list of available taxonomies in this format is available at [Taxonomy_list.md](https://github.com/AllenInstitute/scrattch.taxonomy/blob/main/Taxonomy_list.md)**. As of 10 November 2023, these represent published single-cell RNAseq taxonomies from the Allen Institute and other groups, and are largely complementary to taxonomies included at the [Brain Knowledge Platform](https://portal.brain-map.org/atlases-and-data/bkp), although efforts to integrate are ongoing.

## Documentation

You can find a detail description of all scrattch.taxonomy functions here: ![Documentation](https://github.com/AllenInstitute/scrattch.taxonomy/blob/main/scrattch.taxonomy_0.1.pdf)

Update notes are here: ![Versions](https://github.com/AllenInstitute/scrattch.taxonomy/blob/main/VERSIONS.md)

## Docker

We have setup a docker environemnt for scrattch.taxonomy and scrattch.mapping that contains all the required dependencies and the current version of both scrattch.taxonomy and scrattch.mapping. This docker is accessible through docker hub via: `kapen/scrattch_mapping:0.96.3`.

#### HPC usage:

##### Non-interactive
`singularity shell --cleanenv docker://njjai/scrattch_mapping:0.6.2 Rscript YOUR_CODE.R`

##### Interactive
`singularity shell --cleanenv docker://njjai/scrattch_mapping:0.6.2`


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

2. [**Build and map against a small mouse PatchSeq taxonomy**](https://github.com/AllenInstitute/scrattch.taxonomy/blob/main/examples/build_patchseq_taxonomy.md) This example provides the basics for updating a taxonomy to be compatible with patchseq style mapping and visualization on (internal Allen Institute) MolGen Shiny tools.

3. [**Build and map against a human MTG PatchSeq taxonomy**](https://github.com/AllenInstitute/scrattch.taxonomy/blob/main/examples/build_MTG_patchseq_taxonomy.md) This example provides shows how to create a standard taxonomy and update it to be compatible with patchseq style mapping and visualization on (internal Allen Institute) MolGen Shiny tools. This example essentially combines examples 1 and 2 and applies them to human neocortical data sets.  Data is from Hodge et al. (2019) for human MTG and patch-seq examples are from Berg et al (2021) (layer 2-3, excitatory neurons). 

4. [**Build and map against a MapMyCells taxonomy**](https://github.com/AllenInstitute/scrattch.taxonomy/blob/main/examples/build_taxonomy_MapMyCells.md) This tutorial shows how to run the MapMyCells python mapping algorithm against a taxonomy.
   
## Reporting issues

If you run into any issues, please let Nelson and Jeremy know or [**create a new issue in the 'Issues' tab above**](https://github.com/AllenInstitute/scrattch-taxonomy/issues).

## TODO

- [ ] Update documentation.

## Done
