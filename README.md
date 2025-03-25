# scrattch.taxonomy

Generalized taxonomy building scripts for RNA-seq based taxonomies following the [Allen Institute schema](https://github.com/AllenInstitute/scrattch.taxonomy/tree/main/schema).

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
2. [**Check a Shiny taxonomy and convert to AIT**](https://github.com/AllenInstitute/scrattch.taxonomy/blob/main/examples/checktaxonomy.md) This example shows how you can use `checkTaxonomy` to see what does not abide by the AIT schema, and then use the resulting information to convert to AIT format. *(Note that this example will only work for Allen Institute employees as it points to a private file.)*
3. [**Create a human MTG taxonomy in AIT format with a neuron only 'child' taxonomy**](https://github.com/AllenInstitute/scrattch.taxonomy/blob/main/examples/build_MTG_taxonomy.md) This example provides a step-by-step process for downloading human MTG data from adult neurotypical humans along with the associated SEA-AD taxonomy [(from here)](https://portal.brain-map.org/atlases-and-data/rnaseq/human-mtg-10x_sea-ad), converting it to an AIT file that aligns with the AIT schema, and adding a child taxonomy subsetting to only neuronal types for use with Patch-seq mapping (see scrattch.patchseq library).
   
## Reporting issues

If you run into any issues, please let Nelson and Jeremy know or [**create a new issue in the 'Issues' tab above**](https://github.com/AllenInstitute/scrattch-taxonomy/issues).

## TODO

- [ ] Update documentation.

## Done
