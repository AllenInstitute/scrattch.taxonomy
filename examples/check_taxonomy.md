# Tutorial: Validate a Shiny taxonomy 

#### Required inputs:

#### Additional prerequisites:

#### Check taxonomy:

Use: `singularity shell --cleanenv docker://njjai/scrattch_mapping:0.9.1`

```R
## Load scrattch.taxonomy
library(scrattch.taxonomy)

taxonomy_dir = "/allen/programs/celltypes/workgroups/rnaseqanalysis/HMBA/Aim1_Atlases/Macaque/BasalGanglia/AIT"
taxonomy_file = "Macaque_basalganglia_AIT.h5ad"

ait.anndata = loadTaxonomy(taxonomy_dir, anndata_file=taxonomy_file)

ait.anndata$uns$valid = checkTaxonomy(ait.anndata)
```
