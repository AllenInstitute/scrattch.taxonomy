![image](https://github.com/AllenInstitute/scrattch.taxonomy/assets/25486679/abbd552a-e490-40c0-b665-5d50c960bde0)# AIT / CAS / BKP schema integration

Several competing schema have been created for packaging of taxonomies, data sets, and associated metadata and annotations.  This document aims to align three such schema and propose a way of integrating them into the Allen Institute Taxonomies (AIT) .h5ad file format presented as part of this GitHub repository. The three standards are:

1. **AIT** (described herein)
2. **[Cell Annotation Schema](https://github.com/cellannotation/cell-annotation-schema/) (CAS)**: this schema is becoming more widely-used in the cell typing field as a whole because it is largely compatible with [the CZ CELLxGENE schema](https://github.com/chanzuckerberg/single-cell-curation/blob/main/schema/3.0.0/schema.md). It is also compabible with [Cell Annotation Platform](https://celltype.info/) (CAP) and with [Taxonomy Development Tools](https://brain-bican.github.io/taxonomy-development-tools/) (TDT). CAS has both a general schema and a BICAN-associated schema, both of which are considered herein.
3. **Brain Knowledge Platform (BKP)**: this schema isn't publicly laid out anywhere that I can find, but this is the data model used for [Jupyter Notebooks](https://alleninstitute.github.io/abc_atlas_access/intro.html) associated with the [Allen Brain Cell (ABC) Atlas](https://portal.brain-map.org/atlases-and-data/bkp/abc-atlas).  More generally, any data sets to be included in ABC Atlas, [MapMyCells](https://portal.brain-map.org/atlases-and-data/bkp/mapmycells), or other related BKP resources will eventually need to conform to this format.

It's worth noting that all of these schema are still under development, and we hope they will approach a common schema.

## Cell type taxonomy organization

One major challenge in creating a cell type taxonomy schema is in definition of terms such as "taxonomy," "dataset," "annotation," "metadata," and "data."  Our current working model of a taxonomy is shown below; however, it is becoming increasingly important to (at minimum) separate out the data from the other components, and (ideally) separately out all components to avoid the need to download, open, or upload huge and unweildy files and to integrate with under-development databases. 

![Taxonomy_overview](https://github.com/AllenInstitute/scrattch.taxonomy/assets/25486679/9d36e6bc-db14-4d73-8011-23026756ec08)

That said, it is still important for many use cases to have an option of including all of the information listed above in a single h5ad file for use with CELLxGENE, [scrattch.mapping](https://github.com/AllenInstitute/scrattch.mapping), and other analysis tools, and for ease of sharing in a single file format. To this end, we divide the schema components in this document by the following broad category terms (described below), and for each term indicate which component of the schematic in which it can be found.  In cases where the same information is saved in multiple ways in different schema, we group these fields together as well.

**EXAMPLE**:fire::fire::fire:

We highlight areas of the schema that are not currently in sync and that need continued effort using fire emojis, as shown here.  

### Broad category terms

Here are the current categories that all fields are placed in as a starting point for discussion.  Cases where the location of a field is unclear or ambiguous are called out :fire: .

* **[Data](https://github.com/AllenInstitute/scrattch.taxonomy/blob/KL_div/schema/aligned_schema.md#data)**: This includes anything critical for understanding the cell by gene matrix and to link it with other components.  This includes data (raw and processed), gene information, and cell identifiers.  *Note that for the purposes of this schema, we are excluding raw data (fastq, bam files, etc.) from consideration and are starting from the count matrix.*
* **[Assigned metadata](https://github.com/AllenInstitute/scrattch.taxonomy/blob/KL_div/schema/aligned_schema.md#assigned-metadata)**: This includes cell-level metadata that is assigned at some point in the process between when a cell goes from the donor to a value in the data, and (in theory) can be ENTIRELY captured by values in Allen Institute, BICAN, or related standardized pipelines.  It includes things like donor metadata, experimental protocols, dissection information, RNA QC metrics, and sequencing metadata.
* **[Calculated metadata](https://github.com/AllenInstitute/scrattch.taxonomy/blob/KL_div/schema/aligned_schema.md#calculated-metadata)**: This includes any cell-level or cluster-level metadata that can be calculated explicitly from the **Data** and **Assigned Metadata** without the need for human intervention.  Some examples include # reads detected/cell, # UMI/cell, fraction of cells per cluster derived from each anatomic dissections, expressed neurotransmitter genes (quantitatively defined), average QUANTITATIVE_VALUE (e.g., doublet score) per cluster.
* **[Annotations](https://github.com/AllenInstitute/scrattch.taxonomy/blob/KL_div/schema/aligned_schema.md#annotations)**: This includes any fields related to the annotation of clusters or groups of clusters (collectively called "cell sets").  This includes things like cluster levels, cluster relationships, canonical marker genes, links to existing ontologies (e.g., CL, UBERON) based on judgement calls, expert annotations, and dendrograms.
* **[Analysis](https://github.com/AllenInstitute/scrattch.taxonomy/blob/KL_div/schema/aligned_schema.md#analysis)**: This includes any fields included as the result of or required for specific analysis.  Some examples include latent spaces (e.g., UMAP), cluster level gene summaries (e.g., cluster means, proportions), and variable genes.
* **[Tooling](https://github.com/AllenInstitute/scrattch.taxonomy/blob/KL_div/schema/aligned_schema.md#tooling)**: This includes any fields required for specific tools (e.g., cellxgene, TDT, CAS, CAP) that are not strictly part of the taxonomy and that do not fit in any of the above categories.  This includes things like schema versions and redundent fields from above with different column names.

We expect some of these categories to change but feel this is a good starting point.

### Anndata schematic

Within each broad categorical term, fields are ordered by their location in the anndata object: X, layers, obsm, obs, var, uns.

![Schematic](https://github.com/AllenInstitute/scrattch.taxonomy/blob/main/schema/AIT_anndata_schema.png)

And now let's go on to the schema!

# Proposed integrated schema


## Data

This includes anything critical for understanding the cell by gene matrix and to link it with other components.  This includes data (raw and processed), gene information, and cell identifiers.  **Ideally a schema for this will be defined through other BICAN groups, and can be adopted here.**

#### X 

The `X` component contains logCPM normalized expression data (cell x gene).

#### layers

The `layers` component contains the count matrix (cell x gene).

* `counts` :fire::fire::fire: : The count matrix from which `X` was derived in **AIT**. Called `raw.x` in **CELLxGENE**. We should convert `counts` to `raw.x`.

#### obs

The `obs` component contains cell level metadata.

* `cell_label` :fire::fire::fire: : ID corresponding to each individual cell.  This ID MUST be included in the **data** and in every other location to refer to the data (e.g., metadata and annotations). This is also called `cell_id`, `sample_id`, and `sample_name` in various places.  We should align on a single term. 

#### var

The `var` component contains gene level information.

* `gene` :fire::fire::fire: : A vector of gene symbols.  Broadly useful in the community for defining genes but occasionally problematic; called `gene_symbol` in **BKP**.  More generally, some alignment is needed about whether these or the `ensembl_id` are used for the gene identifier column (CELLxGENE uses a very specific version of `ensembl_id` for the INDEX).
* `ensembl_id` :fire::fire::fire: : A vector of corresponding Ensembl IDs for each gene symbol. This is required for disambiguation of gene symbols; called `gene_identifier` in **BKP**. This is optional for **AIT** (but maybe it sholdn't be?).
* `biotype`: biotype from the gtf file (e.g., protein_coding); used in BKP and BICAN for filtering of genes (but optional elsewhere)
* `name`: Longer gene name from the gtf file; used in BKP
* `[additional gene info]`: Optional uncontrolled gene info; could include gene length, Genecode IDs, NCBI identifiers, etc.

#### uns

The `uns` component contains more general information and fields with formatting incompatible with the above components.

* `dataset_metadata` :fire::fire::fire: : TBD information about the data set itself. A standard on this is not established (as far as I know?) but could include some combination of these pieces of information recorded for **Annotations** below: description, dataset_url, title, dataset_doi, author_list, author_name, author_contact, orcid, etc.


## Assigned metadata

This includes cell-level metadata that is assigned at some point in the process between when a cell goes from the donor to a value in the data, and (in theory) can be ENTIRELY captured by values in Allen Institute, BICAN, or related standardized pipelines.  It includes things like donor metadata, experimental protocols, dissection information, RNA QC metrics, and sequencing metadata.  **Ideally a schema for this will be defined through other BICAN groups, and can be adopted here.**

#### obs

The `obs` component contains cell level metadata from the experiment

* `cell_label`: ID corresponding to each individual cell.  See above.
* `[additional cell ID columns]`: Optional additional IDs per cell.  They are not used for taxonomy efforts.  This could include things like IDs for RNA wells, barcodes, or other tracking IDs used for data processing.
* `feature_matrix_label`: ID of the associated feature matrix where the data is stored (if not included in this file). Used in **BKP** when data is found elsewhere for connected cell to data file. 
* `dataset_label` :fire::fire::fire: : Link between each cell and each dataset in **BKP**.  Need clarification on how this differs from feature_matrix_label; for **CAS** this is a taxonomy-level variable in uns called `dataset_url` (I think).  We should align on this too.
* `[COLUMN_NAME]_color` :fire::fire::fire: : Color vector for metadata/taxonomy values in format [COLUMN_NAME]_label. This is ONLY used for molgen-shiny plots, but because of this, some metadata files come with these and some don't and that could cause challenges.  Should revisit how to store colors and how to deal with metadata in both formats.  Should also agree on a standard for which way is preferred.
* `[COLUMN_NAME]_id` :fire::fire::fire: : Same as above, but in this case for the order of metadata values (e.g., the levels of a factor, or ascending order of a numeric)

The `obs` component also contains **experiment metadata** per cell

* `assay` and `assay_ontology_term_id` :fire::fire::fire: : In **CELLxGENE** these correspond to a human-readable modality along with the associated EFO ontology term. We often use the term `modality` in place of `assay` (e.g., 'Smart-seq2'corresponds to 'EFO:0008931', '10x 3' v3'corresponds to 'EFO:0009922'). This is called "library method" in **BKP**. Ideally we will agree on a term for this and it will be provided upstream from BICAN.
* `suspension_type` :fire::fire::fire: : Either "cell", "nucleus", or "na" in **CELLxGENE**. Called `entity` in the **BKP**. We should pick one to use.
* `[batch_condition_columns]`: Zero or more vectors of metadata associated with batches. These are not required, but called out separately by cellxgene for analysis purposes.
* `[additional uncontrolled metadata]`: Additional uncontrolled cell metadata. These are not required, but any additional columns are allowed by all h5ad formats.

The `obs` component also contains **brain region metadata** per cell, but this is still an active area of development

* `brain_region` :fire::fire::fire: : Brain region(s) sampled. Called `tissue_ontology_term_id` in cellxgene; cell_set structures also defined below; called `region_of_interest_label` and `anatomic_division_label` in BKP. Also associated are acronymns, labels, etc.;  More generally need to arrive at a way of dealing with brain regions.  Note that this slot in the **Assigned metadata** is meant to deal with cell-level assignments for brain region (e.g., dissection) and NOT cell set summarizations by brain region, which are included below.
* `tissue` and `tissue_ontology_term_id` :fire::fire::fire: : Along with "tissue" field, these correspond to UBERON terms for the 'brain region' fields that we have (e.g., 'brain' = 'UBERON_0000955').  In process: we need to discuss how to integrate Allen reference brain atlases for mouse and human.

The `obs` component also contains **donor level metadata** per cell

* `donor_id` :fire::fire::fire: : Identifier for the unique individual, ideal from the specimen portal (or other upstream source). This is called `donor_label` in the **BKP**. Should converge on a standard term. More than one identifier may be needed, but ideally for the analysis only a single one is retained and stored here.
* `species` :fire::fire::fire: : Species sampled. This is split into two fields in CAP/cellxgene/BICAN: `organism` (e.g., homo sapiens) and  `organism_ontology_term_id` (e.g., 'NCBITaxon:10090'). For consistency, we should change `species` to `organism` and could write a function to automatically identify the ontology term (I think [GeneOrthology](https://github.com/AllenInstitute/GeneOrthology/blob/main/README.md) already has one).
* `age` :fire::fire::fire: : Currently a free text field for defining the age of the donor. In **CELLxGENE** this is recorded in `development_stage_ontology_term_id` and is HsapDv if human, MmusDv if mouse.  I'm not sure what this means, but more generally, we should align with BICAN on how to deal with this value.
* `sex` :fire::fire::fire: : Placeholder for donor sex. Called `sex_ontology_term_id` (e.g., PATO:0000384/383 for male/female) in **CELLxGENE** and called "donor_sex" in BKP. We should align on a single term.
* `donor_genotype`: One (or sometimes more) column related to the genotype of the animal (for transgenic mice, in particular). Not used for humans and most NHP.
* `self_reported_ethnicity_ontology_term_id`: Controversial field that is required for **CELLxGENE** but otherwise not used. HANCESTRO term if human and 'na' if non-human.
* `disease` and `disease_ontology_term_id`: A human-readable name for a disease and the associated MONDO ontology term (or PATO:0000461 for 'normal'). Used in **CELLxGENE** and ideally we can also adopt for **SEA-AD** and other use cases.

#### uns

The `uns` component contains more general information and fields with formatting incompatible with the above components.

* `assigned_metadata_metadata` :fire::fire::fire: : TBD information about the assigned_metadata itself. This likely is not needed or should be renamed.
* `batch_condition`: List of obs fields that define “batches”; Used by **CELLxGENE** if provided, but otherwise not needed.


## Calculated metadata

This includes any cell-level or cluster-level metadata that can be calculated explicitly from the **Data** and **Assigned Metadata** without the need for human intervention.  Some examples include # reads detected/cell, # UMI/cell, fraction of cells per cluster derived from each anatomic dissections, expressed neurotransmitter genes (quantitatively defined), average QUANTITATIVE_VALUE (e.g., doublet score) per cluster. **Currently none of these are required for the schema, but they are sometimes used for annotation.**

#### obs

The `obs` component contains cell level metadata from the experiment

* `cell_label`: ID corresponding to each individual cell.  See above.
* `[additional uncontrolled metadata]`: Additional uncontrolled cell metadata. These are not required, but any additional columns are allowed by all h5ad formats.
* `feature_matrix_label`, `dataset_label`, `[COLUMN_NAME]_color`, `[COLUMN_NAME]_id`: See above.

#### uns

The `uns` component contains more general information and fields with formatting incompatible with the above components.

* `calculated_metadata_metadata` :fire::fire::fire: : TBD information about the calculated_metadata itself. This likely is not needed or should be renamed.


## Annotations

This includes any fields related to the annotation of clusters or groups of clusters (collectively called "**cell sets**").  This includes things like cluster levels, cluster relationships, canonical marker genes, links to existing ontologies (e.g., CL, UBERON) based on judgement calls, expert annotations, and dendrograms.

#### obs

The `obs` component contains cell level metadata, as above.

* `cell_label`: ID corresponding to each individual cell.  See above.
* `cluster` :fire::fire::fire: : This is the **CRITICAL** column used for cluster annotations. It is the baseline for the majority of cell_annotation columns discussed below. Sometimes called `cluster_label`; There is also an additional `cluster_alias` column used in mouse whole brain data and for BKP that I'm not sure how to wrap in. It's also used for cirrocumulus.  Some discussion might be needed on whether this is one or more columns, and which one is the source of truth. It's also worth noting that this is a prerequisite for annotations, so maybe it better fits in a different category (analysis?).
* `[additional uncontrolled metadata]`: Additional uncontrolled cell metadata. These are not required, but any additional columns are allowed by all h5ad formats.
* `feature_matrix_label`, `dataset_label`, `[COLUMN_NAME]_color`, `[COLUMN_NAME]_id`: See above.

The `obs` component also contains **cluster-level metadata** summarize at the cell level. In theory a majority of this information could be stored in the uns (or in separate json files), but for now we have it listed here for consistency with **Cell annotation schema**

* `[cellannotation_set]` :fire::fire::fire: : This is where **taxonomy levels** are stored in **CAS**. The column name is a string (e.g., "subclass") and the values are `cell_labels` (e.g., "SST").  The equivalent in **BKP** are `Cluster Annotation Term Set`s and in **CAP** is `label_sets`.  This also encapsulates the concept of `cell_ids` from CAP/CAS, since in the h5ad file each row corresponds to a cell and therefore you get the `cell_label` --> `cell_id`s mapping for free.  This is stored as separate files in both **BKP** and **TDT** and so we should confirm appropriate translations and agree on terms for this.  Finally, this concept is critical for scrattch.taxonomy and scrattch.mapping to work properly, but none of the `[cellannotation_set]--XXXX` fields below are needed for **AIT** (although keeping this is best practice!).
* `[cellannotation_set]--cell_set_accession` :fire::fire::fire: : ID corresponding to the cell_set; called the "Cluster Annotation Term" in **BKP**.  Some work still needed on deciding what to name this (CCNXXXXX?, a hash code?, something else?) and HOW to name this (automatically? if so, but what authority).
* `[cellannotation_set]--cell_set_label` :fire::fire::fire: : 
* `[cellannotation_set]--cell_fullname` :fire::fire::fire: : 
* `[cellannotation_set]--cell_ontology_exists` :fire::fire::fire: : 
* `[cellannotation_set]--cell_ontology_term_id` :fire::fire::fire: : 
* `[cellannotation_set]--cell_ontology_term` :fire::fire::fire: : 
* `[cellannotation_set]--rationale` :fire::fire::fire: : 
* `[cellannotation_set]--rationale_dois` :fire::fire::fire: : 
* `[cellannotation_set]--marker_gene_evidence` :fire::fire::fire: : 
* `[cellannotation_set]--canonical_marker_genes` :fire::fire::fire: : 
* `[cellannotation_set]--synonyms` :fire::fire::fire: : 
* `[cellannotation_set]--category_XXXX` :fire::fire::fire: : 
* `[cellannotation_set]--parent_cell_set_accession` :fire::fire::fire: : 

(THIS STUFF GOES ABOVE.  I JUST DON'T WANT TO LOSE IT!)
Cluster label or list of underlying cluster labels (deprecated)	Used to be required for CCN, but can probably be dropped
Longer name (e.g., "Somatostatin interneuron")	CCN: either this or cell_label is "cell_set_preferred_alias"
True/false call about whether a cell ontology term exists	(This seems redundant to me with next two rows)
Highest resolution Cell Ontology term (ID)	CCN: "cell_set_ontology_tag"
Highest resolution Cell Ontology term (name)	CCN: "cell_set_structure"/"cell_set_aligned_alias"
Free text evidence for cell annotations	
Comma-separated publication DOI's of rationale	CCN: "cell_set_alias_citation"
Comma-separated marker genes used as evidence for cell type annotation (e.g., by NS-Forest)	This is reserved for ontology markers.  See var below for general marker genes.
Comma separated list of canonical marker genes	I don't understand how this differs from "marker_gene_evidence"
Comma-separated aliases (e.g., "neuroglial cell, glial cell, neuroglia")	CCN: "cell_set_additional_alias"
[Same as cell_XXXX above for several fields, e.g., fullname, cell_ontology_term, etc.]	(I'm not sure how this differs from the cell_XXXX)
ID corresponding to the parent cell_set	ID for each parent cell set; possibly the parent "Cluster Annotation Term" in knowlegebase?




## Analysis

#### obsm

The `obsm` component contains all dimensionality reductions of the taxonomy (cell x dim). For all fields listed below, columns are of the format '[FIELD]_#' where # is 1, 2, 3, etc..

* `umap` :fire::fire::fire: : 2 (or more)-dimensional representation of cells in **AIT**. Must be of the form `X_[...]` for use with **CELLxGENE**.  Only the first two dimensions are used for AIT and CELLxGENE, but 3 dimensions can be used for cirrocumulus.
* `pca`: Additional terms for embedding multi-dimensional principal components and latent spaces
* `scVI`: Additional terms for embedding multi-dimensional principal components and latent spaces

#### obs

The `obs` component contains cell level metadata including cell type annotation, brain region, species and additional taxonomy specific fields. Standardized nomenclature for cell metadata to be decided.

* `cell_label`: ID corresponding to each individual cell.  See above.
* `brain_region`
* `species`
* `age`
* ...

#### var

The `var` component contains gene level metadata including which genes are highly variable or marker genes. Additional information may enclude Ensemble IDs.

* `gene`: A vector of gene symbols
* `ensembl_id`: A vector of corresponding Ensembl IDs for each gene symbol. (**Optional**)
* `highly_variable_genes`: A logical vector (T/F) indicating which genes are highly variable.
* `marker_genes`: A logical vector (T/F) indicating which genes are markers used to build dendrogram. (**Optional**)


#### uns

The `uns` component contains taxonomy associated files useful for reproducing analysis or mapping against the taxonomy.

* `dend`: A file path to the dendrogram .RData object. Will be replaced by language agnostic encoding of dendrogram.
* `filter`: A logical vector indicating which cells to filter out for a specific taxonomy version.
* `QC_markers`: Marker gene expression in on-target and off-target cell populations, useful for patchseq analysis.
* `mode`: Taxonomy mode that determins which `filter` to use. See scrattch.mapping documentation.
* `clustersUse`: A vector of cluster names to use for taxonomy.
* `clusterInfo`: A data.frame of cluster information.
* `taxonomyName`: Taxonomy name and AIT identifier.
* `taxonomyDir`: Directory where taxonomy related files can be located. (**Optional**)
