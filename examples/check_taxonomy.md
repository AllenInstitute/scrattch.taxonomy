# Tutorial: Validate a Shiny taxonomy 

#### Required inputs:

#### Additional prerequisites:

### Check taxonomy:

Use: `singularity shell --cleanenv docker://njjai/scrattch_mapping:1.0.0`

```R
## Load scrattch.taxonomy
library(scrattch.taxonomy)

taxonomy_dir = "/allen/programs/celltypes/workgroups/rnaseqanalysis/HMBA/Aim1_Atlases/Macaque/BasalGanglia/AIT"
taxonomy_file = "Macaque_basalganglia_AIT.h5ad"

ait.anndata = loadTaxonomy(taxonomy_dir, anndata_file=taxonomy_file)

ait.anndata = checkTaxonomy(ait.anndata, print.messages=TRUE)
```

This is an old file which contains most of the schema elements, but is missing a bunch of other elements.  **If you downloaded the file yourself or created it with docker versions 0.9 or higher, your file should pass and you can probably stop here.**

---

### Convert file to appropriate AIT format

This section steps through the process of building an updated AIT-format file using the data in ait.anndata above. 

First, let's walk through each of the warnings and errors presented above.

```R
taxonomy.anno <- ait.anndata$obs

# WARNING: AIT.anndata$X has high values.  Please confirm this is log-normalized.
# This error will be resolved with buildTaxonomy
taxonomy.counts <- t(ait.anndata$X)

# WARNING: the following AIT.anndata$uns columns are **REQUIRED** for the schema: hierarchy, filter, cluster_info, schema_version.
hierarchy = list("cluster_id") # This taxonomy does not have a hierarchy yet, but we can still put in level 1.
# filter, cluster_info, and schema_version are generated automatically

# WARNING: AIT.anndata$var does not contain highly_variable_genes[_name], which is recommended for generating UMAPs and dendrograms.
# WARNING: AIT.anndata$var does not contain marker_genes[_name], which is recommended for generating UMAPs and dendrograms.
binary.genes = top_binary_genes(taxonomy.counts, taxonomy.anno$cluster_id, 1000)

# The anndata.obs element: ensembl_id contains 0% ensembl terms.
# Common NCBI taxIDs: Human = 9606; Mouse = 10090; Macaque (rhesus) = 9544; Marmoset = 9483
ensembl_id <- geneSymbolToEnsembl(gene.symbols = rownames(ait.anndata$var), ncbi.taxid = 9606)

# WARNING: AIT.anndata$obsm does not include any embeddings, which are required for some downstream functions.
pcs  <- prcomp(logCPM(taxonomy.counts)[binary.genes,], scale = TRUE)$rotation
umap.coords = umap(pcs[,1:30])$layout
rownames(umap.coords) = colnames(taxonomy.counts)


# ERROR: taxonomy modes with filters are not found. Allen Institute Taxonomy requires at least a standard mode with all cells included which should have been created with `buildTaxonomy`.  Likely this h5ad is an earlier version of Allen Institute Taxonomy (AIT) format and should be remade.
# This error will be resolved with buildTaxonomy
```

Now lets check (and if needed update) the taxonomy metadata. This will look for common misspellings in schema columns and check or add relevant ontology terms for "organism", "anatomical_region", "self_reported_sex", "self_reported_ethnicity", "assay", and "disease". **Note that these functions will CHANGE your metadata file! Please review the results carefully before proceeding.**

```R
# First do some basic QC of the metadata
full.taxonomy.anno <- updateTaxonomyMetadata(taxonomy.anno, print.messages=TRUE)

# Now, compute ontology terms for relevant columns
full.taxonomy.anno <- computeOntologyTerms(full.taxonomy.anno$metadata, print.messages=TRUE)

# Clearly some of these changes are correct and some are incorrect. For now, we will incorporate the subset of columns that look correct.
#taxonomy.anno <- full.taxonomy.anno

update_CN <- c("anatomical_region_ontology_term_id", 
               "brain_region_ontology_term_id", 
               "self_reported_sex",
               "self_reported_sex_ontology_term_id", 
               "assay",
               "assay_ontology_term_id")

for(cn in update_CN) taxonomy.anno[,cn] <- full.taxonomy.anno$metadata[,cn]

rownames(taxonomy.anno) = colnames(taxonomy.counts)
```

Finally, regenerate the AIT-formatted file using the above information

```R
## Build Allen Insitute Taxonomy, for large taxonomies you can pass in tpm and cluster_stats if pre-computed.
AIT.anndata = buildTaxonomy(meta.data = taxonomy.anno,
                            title = "AIT.anndata",
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
                            subsample=100,
                            ##
                            add.dendrogram.markers = FALSE, # This fails for some reason on this taxonomy, but normally this is set to TRUE
                            addMapMyCells = FALSE           # This fails for some reason on this taxonomy, but normally this is set to TRUE
                            )

## Check whether the UPDATED taxonomy file is valid
AIT.anndata$uns$valid = checkTaxonomy(AIT.anndata)
```