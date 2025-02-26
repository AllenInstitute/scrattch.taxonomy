# scrattch.taxonomy

Generalized taxonomy building scripts for RNA-seq based taxonomies following the [Allen Institute schema]([https://github.com/AllenInstitute/scrattch.taxonomy/tree/main/schema).

**A list of available taxonomies in this format is available at [Taxonomy_list.md](https://github.com/AllenInstitute/scrattch.taxonomy/blob/main/Taxonomy_list.md)**. As of 10 November 2023, these represent published single-cell RNAseq taxonomies from the Allen Institute and other groups, and are largely complementary to taxonomies included at the [Brain Knowledge Platform](https://portal.brain-map.org/atlases-and-data/bkp), although efforts to integrate are ongoing.

## Documentation

You can find a detail description of scrattch.taxonomy (OUTDATED -- FIX) functions here: ![Documentation](https://github.com/AllenInstitute/scrattch.taxonomy/blob/main/scrattch.taxonomy_0.1.pdf)

Update notes are here: ![Versions](https://github.com/AllenInstitute/scrattch.taxonomy/blob/main/VERSIONS.md)

## Installation

### Using docker (recommended)
We have setup a docker environemnt for scrattch.taxonomy, scrattch.mapping, and scrattch.patchseq that contains all the required dependencies and the current version of all scrattch packages. **See [the readme](https://github.com/AllenInstitute/scrattch/blob/master/README.md#using-docker) for [the parent scrattch package](https://github.com/AllenInstitute/scrattch) for the most up-to-date docker information.**

### Directly from GitHub (strongly discouraged)

While we advise using the provided docker, you can also install scrattch.taxonomy directly from GitHub as follows:

```
devtools::install_github("AllenInstitute/scrattch.taxonomy")
```

This strategy **might not work** due to complicated dependencies. Also note that `doMC` may need to be installed manually from [HERE](https://r-forge.r-project.org/R/?group_id=947) if you use Windows. Vignettes are provided below.

## Usage examples

1. [**Build a Shiny taxonomy**](https://github.com/AllenInstitute/scrattch.taxonomy/blob/main/examples/build_taxonomy.md) This example provides the basics for creating a new taxonomy compatible with scrattch.mapping mapping functions and (internal Allen Institute) MolGen Shiny tools.
   
## Reporting issues

If you run into any issues, please let Nelson and Jeremy know or [**create a new issue in the 'Issues' tab above**](https://github.com/AllenInstitute/scrattch-taxonomy/issues).

## TODO

- [ ] Update documentation.

## Done
