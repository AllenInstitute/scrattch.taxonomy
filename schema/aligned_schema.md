# AIT / CAS / BKP schema integration

*(Note: An evolving version of this standard is available **[as a Google Doc](https://docs.google.com/document/d/1nj6LHUPoo3JnNwZ7PTdniT9pBPsoJr1B/edit?usp=sharing&ouid=113573359044104089630&rtpof=true&sd=true)**).*

Several competing schema have been created for packaging of taxonomies, data sets, and associated metadata and annotations.  This document aims to align three such schema and propose a way of integrating them into the Allen Institute Taxonomies (AIT) .h5ad file format presented as part of this GitHub repository. The three standards are:

1. **AIT** (described herein)
2. **[Cell Annotation Schema](https://github.com/cellannotation/cell-annotation-schema/) (CAS)**: this schema is becoming more widely-used in the cell typing field as a whole because it is largely compatible with [the CZ CELLxGENE schema](https://github.com/chanzuckerberg/single-cell-curation/blob/main/schema/3.0.0/schema.md). It is also compabible with [Cell Annotation Platform](https://celltype.info/) (CAP) and with [Taxonomy Development Tools](https://brain-bican.github.io/taxonomy-development-tools/) (TDT). CAS has both a general schema and a BICAN-associated schema, both of which are considered herein.  CAS can be embedded in the header (`uns`) of an AIT/Scraatch.taxonomy file, where it functions as a store of extended information about an annotation, including ontology term mappings, evidence for annotation (from annotation transfer and marker expression).
3. **Brain Knowledge Platform (BKP)**: this schema isn't publicly laid out anywhere that I can find, but this is the data model used for [Jupyter Notebooks](https://alleninstitute.github.io/abc_atlas_access/intro.html) associated with the [Allen Brain Cell (ABC) Atlas](https://portal.brain-map.org/atlases-and-data/bkp/abc-atlas).  More generally, any data sets to be included in ABC Atlas, [MapMyCells](https://portal.brain-map.org/atlases-and-data/bkp/mapmycells), or other related BKP resources will eventually need to conform to this format.

It's worth noting that all of these schema are still under development, and we hope they will approach a common schema.

[Taxonomy_field_mappings](https://docs.google.com/spreadsheets/d/1PhsOipO0yCrtTGrkWXLU2Tj2m4qFPgdKID0SqetYN0Y/edit#gid=0) in table form.

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

Here is a graphical representation of these terms in the context of data, metadata, and taxonomies:
![image](https://github.com/AllenInstitute/scrattch.taxonomy/assets/25486679/eaf6b3d3-0b5f-49fc-9a49-2b7168605964)


### Anndata schematic

Within each broad categorical term, fields are ordered by their location in the anndata object: X, layers, obsm, obs, var, uns.

![Schematic](https://github.com/AllenInstitute/scrattch.taxonomy/blob/main/schema/AIT_anndata_schema.png)

And now let's go on to the schema!

# Proposed integrated schema


## Data

This includes anything critical for understanding the cell by gene matrix and to link it with other components.  This includes data (raw and processed), gene information, and cell identifiers.  **Ideally a schema for this will be defined through other BICAN groups, and can be adopted here.**

#### X 

The `X` component contains logCPM normalized expression data (cell x gene).

#### layers → (Use raw.X instead)

The `layers` component contains the count matrix (cell x gene).

* `counts` :fire::fire::fire: : The count matrix from which `X` was derived in **AIT**. Called `raw.x` in **CELLxGENE**. We should convert `counts` to `raw.x`.

#### obs

The `obs` component contains cell level metadata.

* `cell_label` :fire::fire::fire: : **Use `cell_id`** ID corresponding to each individual cell.  This ID MUST be included in the **data** and in every other location to refer to the data (e.g., metadata and annotations). This is also called `cell_id`, `sample_id`, and `sample_name` in various places.  We should align on a single term. 

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

* `assay` and `assay_ontology_term_id` :fire::fire::fire: : In **CELLxGENE** these correspond to a human-readable modality along with the associated EFO ontology term. We often use the term `modality` in place of `assay` (e.g., 'Smart-seq2'corresponds to 'EFO:0008931', '10x 3' v3'corresponds to 'EFO:0009922'). This is called "library method" in **BKP**. Ideally we will agree on a term for this and it will be provided upstream from BICAN. Called `Modality` in taxonomy Google Sheet.
* `suspension_type` :fire::fire::fire: : Either "cell", "nucleus", or "na" in **CELLxGENE**. Called `entity` in the **BKP**. We should pick one to use.
* `[batch_condition_columns]`: Zero or more vectors of metadata associated with batches. These are not required, but called out separately by cellxgene for analysis purposes.
* `[additional uncontrolled metadata]`: Additional uncontrolled cell metadata. These are not required, but any additional columns are allowed by all h5ad formats.

The `obs` component also contains **brain region metadata** per cell, but this is still an active area of development

* `brain_region` :fire::fire::fire: : Brain region(s) sampled. Called `tissue_ontology_term_id` in cellxgene; cell_set structures also defined below; called `region_of_interest_label` and `anatomic_division_label` in BKP. Also associated are acronymns, labels, etc.;  More generally need to arrive at a way of dealing with brain regions.  Note that this slot in the **Assigned metadata** is meant to deal with cell-level assignments for brain region (e.g., dissection) and NOT cell set summarizations by brain region, which are included below.
* `tissue` and `tissue_ontology_term_id` :fire::fire::fire: : Along with "tissue" field, these correspond to UBERON terms for the 'brain region' fields that we have (e.g., 'brain' = 'UBERON_0000955').  In process: we need to discuss how to integrate Allen reference brain atlases for mouse and human.

The `obs` component also contains **donor level metadata** per cell

* `donor_id` :fire::fire::fire: : Identifier for the unique individual, ideal from the specimen portal (or other upstream source). This is called `donor_label` in the **BKP**. Should converge on a standard term. More than one identifier may be needed, but ideally for the analysis only a single one is retained and stored here.
* `species` :fire::fire::fire: : Species sampled. This is split into two fields in CAP/cellxgene/BICAN: `organism` (e.g., homo sapiens) and  `organism_ontology_term_id` (e.g., 'NCBITaxon:10090'). For consistency, we should change `species` to `organism` and could write a function to automatically identify the ontology term (I think [GeneOrthology](https://github.com/AllenInstitute/GeneOrthology/blob/main/README.md) already has one). Called `Species name` and `Species ID` in taxonomy Google Sheet. 
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
*  `cell_annotation_schema` - extended metadata about annotations and labelsets stores in JSON.
   * extended metadata about labelsets is stored under the `labelsets' key` in CAS (master documentation in [CAS - BICAN extension](https://github.com/cellannotation/cell-annotation-schema/blob/main/build/BICAN_schema.md) under `labelsets`.  Summarised here:
      * name (string, required): name of annotation key. This corresponds to the `obs` key name e.g. SubClass, Neurotransmitter.  In `BKP` this corresonds to cluster_annotation_term_set.name. 
      * description (string): Some text describing what types of cell annotation this annotation key is used to record, e.g. This labelset is used to record neurotransmitter.
      * annotation_method (string): The method used for creating the cell annotations. This MUST be one of the following strings: 'algorithmic', 'manual', or 'both' . Must be one of: ["algorithmic", "manual", "both"].
      * annotation_rank: An integer used to indicate hierarchy level with 0 being the most granular.  e.g. a hierarchy consisting of cluster, subtype, type, class would have ranks 0,1,2,3 respectively.  For non hierarchical labelsets, this should be left blank.  BKP level is similar but reversed (leaf nodes are the most granular).
   * extended metadata about cell sets is stored under the `annotations' key in CAS. (master documentation in [CAS - BICAN extension](https://github.com/cellannotation/cell-annotation-schema/blob/main/build/BICAN_schema.md) under `labelsets`.  The values of labelset and cell_label corresponds to obs key value pairs:
   * labelset: This corresponds to the `obs` key name e.g. SubClass, Neurotransmitter.  The equivalent in **BKP** is `Cluster Annotation Term Set Name' (BKP also has asn ID for this,  referred to as `Cluster Annotation Term Set Label`
   * cell_label: This corresponds to a value of the obs key: e.g. MSN-D1, cholinergic
   * cell_fullname`: The longer name for a cell type, with abbreviations spelled out (e.g., "Somatostatin interneuron 1" rather than "SST 1").  This was called the `cell_set_preferred_alias` in CCN.
   * cell_ontology_term_id`: Highest resolution Cell Ontology term (ID); was called `cell_set_ontology_tag` in CCN
   * cell_ontology_term`: Highest resolution Cell Ontology term (name); was called `cell_set_structure` in CCN and was also largely mapping to the `cell_set_aligned_alias`  
   * `cell_set_accession` :fire::fire::fire: : ID corresponding to the cell_set; called the "Cluster Annotation Term Label" in **BKP**.  ***Some work still needed on deciding what to name this (CCNXXXXX?, a hash code?, something else?) and HOW to name this (automatically? if so, but what authority).***
   * `rationale`: Free text field descrbing the evidence for cell annotations.
   * `rationale_dois` :fire::fire::fire: : List of DOIs of supporting papers.  CAS stores this as a list, but there is some variation in the delimiter used in string representation for reporting and editing purposes (`,`,`|`,`/#/` purposes represented as strings 
    * `marker_gene_evidence`: List of genes used as evidence that this labelset corresponds to the cell type referred to by the cell_label.  For example, expression of DRD1 might be used as evidence that that this cell set correpsonds to DRD1 cells.
    * `synonyms`: Lst of commonly used synonyms/aliases for this cell type (e.g., "neuroglial cell, glial cell, neuroglia"); was called `cell_set_additional_alias` in CCN.
    * `parent_cell_set_accession` :fire::fire::fire: : ID corresponding to the parent cell_set. In BKP this corresponds to  `cluster_annotation_term.parent_term_label` 
    * author_annotation_fields:  CAS supports any additional fields that taxonomy authors may need to associate with cell sets.
    *  `transferred_annotations` :fire::fire::fire: : This key stores a sub-table of annotation transfers: 
       * transferred_cell_label (string): Transferred cell label.
       * source_taxonomy (string): PURL of source taxonomy.
       * source_node_accession (string): accession of node that label was transferred from.
       * algorithm_name (string): 
       * comment (string): Free text comment on annotation transfer.


    




## Annotations

This includes any fields related to the annotation of clusters or groups of clusters (collectively called "**cell sets**").  This includes things like cluster levels, cluster relationships, canonical marker genes, links to existing ontologies (e.g., CL, UBERON) based on judgement calls, expert annotations, and dendrograms.

#### obs

The `obs` component contains cell level metadata, as above.

* In AnnData files, the ID corresponding to each individual cell is stored in the obs index.
* `cluster` :fire::fire::fire: : This is the **CRITICAL** column used for cluster annotations. It is the baseline for the majority of cell_annotation columns discussed below. Sometimes called `cluster_label`; There is also an additional `cluster_alias` column used in mouse whole brain data and for BKP that I'm not sure how to wrap in. It's also used for cirrocumulus.  Some discussion might be needed on whether this is one or more columns, and which one is the source of truth. It's also worth noting that this is a prerequisite for annotations, so maybe it better fits in a different category (analysis?).
* `[additional uncontrolled metadata]`: Additional uncontrolled cell metadata. These are not required, but any additional columns are allowed by all h5ad formats.
* `feature_matrix_label`, `dataset_label`, `[COLUMN_NAME]_color`, `[COLUMN_NAME]_id`: See above.

The `obs` component also contains **cell-level metadata** summarize at the cell level. In theory a majority of this information could be stored in the uns (or in separate json files), but for now we have it listed here for consistency with **Cell annotation schema**



#### var

The `var` component contains gene level metadata.

* `gene`: Same vector included in "data" to link between files.
* `marker_genes_[...]` :fire::fire::fire: : A set of logical vectors (T/F) indicating which genes are markers used to build dendrogram, or for other purposes. The `[...]` part of the name links to additional metadata in the `uns`. This needs to be **UPDATED** in **AIT** to allow multiple marker gene sets; markers currently stored differently in CAP.

#### uns

The `uns` component contains more general information and fields with formatting incompatible with the above components.

* `dend` :fire::fire::fire: : A json formatted dendrogram used for tree mapping. Created by **scrattch.taxonomy** if not provided. Sometimes used for taxonomy annotation, but we are moving away from it with larger taxonomies and so this may now make more sense in the "analysis" category.
* `labelsets` :fire::fire::fire: : ***CRITICAL extra component***; Equilalent to `Cluster annotation term set` in **BKP**. This is saved as a data frame representation (or is a list of data frames needed?), with some information about each `[cellannotation_set]` set of columns (e.g., subclass, class, neurotransmitter, etc.). Specifically: "name", "description", and "rank" (0 most specific) and some information about provenance are needed for each labelset.
* `cell_set_relationships` :fire::fire::fire: : **NEW** proposed mechanism for dealing with sibling relationships for things like gradients, trajectories, constellation diagrams, etc.. This is stored as a data frame (table) of all relations with five columns: cells_set_accession1, cell_set_accession2, relation_label, value, direction.  Could alternatively be stored as a JSON representation that unpacks into a dataframe.
* `filter` :fire::fire::fire: : Indicator of which cells to use for a given child taxonomy (subset), saved as a list of vectors. Each entree in this list is named for the relevant "mode" and has TRUE/FALSE calls indicating whether a cell is filtered out (e.g., the "standard" taxonony is all FALSE). This is critical for how child taxonomies are defined and implemented in **scrattch.taxonomy** but differs from how taxonomies are stored in all other schemas--*some discussion may be needed*.  
* `taxonomyName` :fire::fire::fire: : Taxonomy name (e.g., "AIT30"); called `title` in **cellxgene**, not sure about other schema. Called `Taxonomy short name` in taxonomy Google Sheet. 
* `taxonomy_id` :fire::fire::fire: : Taxonomy ID in CCN format (e.g., "CCN030420240"); TBD how this is generated, but MUST be globally unique. Also used as part of PURL (I think).  Called `Taxonomy ID` in taxonomy Google Sheet too.
* `description` :fire::fire::fire: : Free text description of the taxonomy (or of the dataset on **CAP**). This is also something we are adding as a requirement for the **BKP**, and I think should be required for all taxonomies.  Called `Description` in taxonomy Google Sheet.
* `taxonomy_citation`: "|"-separated publication DOI's of the taxonomy (e.g., "doi:10.1038/s41586-018-0654-5"). Called `Publication` in taxonomy Google Sheet.
* `marker_gene_metadata` :fire::fire::fire: : Data frame of Marker genes x dims that includes metadata for marker gene sets in `var` above; **NEW** and required if `marker_genes_[…]` is provided.  At minimum a name (matching above) and description are needed, but potentially other things (e.g., what is it for, with controlled vocabulary).
* `transferred_annotations_metadata` :fire::fire::fire: : Data frame of info about each transferred annotation column: source_taxonomy, algorithm_name, comment; Still some work on the best way to code this, but it is important.  Linked to data in `var` above.  This is for taxonomy-level metadata. This is also already encoded in **TDT**--how? 
* `taxonomyDir` :fire::fire::fire: : Location of the h5ad file; we might be able to remove this, since it is redundant with `dataset_url` and/or `matrix_file_id`.  Called `Taxonomy file location` in taxonomy Google Sheet.
* `dataset_url` :fire::fire::fire: : PURL of taxonomy; Possibly a redundant field, but critical; also `publication_url` and `cellannotation_url` (unclear how different)
* `matrix_file_id` :fire::fire::fire: : Like `dataset_url`;  e.g. CellXGene_dataset:8e10f1c4-8e98-41e5-b65f-8cd89a887122; Note: needs to be extended to allow for more than one file and connected to `feature_matrix_label` in `obs`.  We need this field!
* `author_list`: List of all collaborators, comma separated [First] [Last]; Useful in general, even though currently only required by **CAP**. Called `Taxonomy Users` in taxonomy Google Sheet.
* `author_name`: The primary author [First Name] [Last Name]; in **CCN** was called "taxonomy_author"; In CCN also seperated by cell_set with "cell_set_alias_assignee"; Called `Point person name` in taxonomy Google Sheet.
* `author_contact`: Valid email address; Called `Point person email` in taxonomy Google Sheet.
* `orcid`: Valid ORCID; Called `Point person ORCID` in taxonomy Google Sheet.
* `annotation_source` : Additional metadata about annotation algorithm; Similar to taxonomy algorithm info stored for CCN


## Analysis

This includes any fields included as the result of or required for specific analysis.  Some examples include latent spaces (e.g., UMAP), cluster level gene summaries (e.g., cluster means, proportions), and variable genes.  These may not need to match between schemas (or even be encoded into schemas).


#### obsm

The `obsm` component contains all dimensionality reductions of the taxonomy (cell x dim). For all fields listed below, columns are of the format '[FIELD]_#' where # is 1, 2, 3, etc..

* `umap` :fire::fire::fire: : 2 (or more)-dimensional representation of cells in **AIT**. Must be of the form `X_[...]` for use with **CELLxGENE**.  Only the first two dimensions are used for AIT and CELLxGENE, but 3 dimensions can be used for cirrocumulus.
* `pca`: Additional terms for embedding multi-dimensional principal components and latent spaces
* `scVI`: Additional terms for embedding multi-dimensional principal components and latent spaces

#### var

The `var` component contains gene level metadata.

* `gene`: Same vector included in "data" to link between files.
* `highly_variable_genes`: A logical vector (T/F) indicating which genes are highly variable. Used for correlation-based mapping in **scrattch.mapping**.
* `marker_genes_[...]`: Potentially additional sets of logical vectors for marker genes, as defined above. 

#### uns

The `uns` component contains taxonomy associated files useful for reproducing analysis or mapping against the taxonomy.

* `dend`: See above. This may fit better here.
* `QC_markers`: Marker gene expression in on-target and off-target cell populations, useful for patchseq analysis.  Also includes information about KL divergence calculations and associated QC calls. Defined by buildPatchseqTaxonomy.
* `filter`: Indicator of which cells to use for a given child taxonomy (subset), as defined above.  
* `mode`: Taxonomy mode that determines which `filter` to use (e.g., that indicates which child taxonomy to map against). Several of the other analysis components of the `uns` have things saved with mode as the name in the h5ad file. See **scrattch.mapping** documentation. Mode is the `Taxonomy short name` in taxonomy Google Sheet for a child taxonomy with the `Parent taxonomy` listed as the `taxonomyName`.
* `clustersUse` :fire::fire::fire: : A vector of cluster names to use for taxonomy. We should be able to remove this 
* `clusterInfo` :fire::fire::fire: : A data.frame of cluster information. We should be able to remove this 
* `marker_gene_metadata`: Metadata about any new marker gene lists added, if any. See above.
* `development_date`: Data of taxonomy development.  Required for Google Sheet.  Potentially not needed if we want to infer from `taxonomy_id`.
* `public`: logical flag indicating whether taxonomy should be public or private. Required for Google Sheet.  Potentially not needed if we want to infer from PURL/GitHub somehow.
* `annotation_sheet`: Link to annotation sheet (ideally a TDT GitHub repo link for communinal annotation).  An optional slot in the Google sheet. I'm not sure if this is listed above somewhere.
* `purpose`: Controlled vocabulary (currently "General" and/or "Patch-seq"). Required for Google Sheet at the moment.


## Tooling

This includes any fields required for specific tools (e.g., cellxgene, TDT, CAS, CAP) that are not strictly part of the taxonomy and that do not fit in any of the above categories.  This includes things like schema versions and redundent fields from above with different column names.  These may not need to match between schemas (or even be encoded into schemas).  We may want to merge this category with Analysis :fire::fire::fire: .

#### obs

The `obs` component contains cell level metadata, as above.

* `cell_label`: ID corresponding to each individual cell.  See above.
* `[cellannotation_set]--parent_cell_set_accession`: ID corresponding to the parent cell_set. If not needed for annotations, definitely needed for tooling.

#### uns

The `uns` component contains taxonomy associated files useful for reproducing analysis or mapping against the taxonomy.

* `schema_version`: cellxgene schema version (e.g., "3.0.0")
* `[...]_color` :fire::fire::fire: : RGB color vector for metadata `[...]`; required only for selecting colors in cirrocumulus.  This may be the same as the `[COLUMN_NAME]_color` column above.
* `cellannotation_schema_version`: **CAS** schema version '[MAJOR].[MINOR].[PATCH]'
* `cellannotation_timestamp` :fire::fire::fire: : Timestamp when published: %yyyy-%mm-%dd %hh:%mm:%ss; Useful in general, even though currently only required by CAP; also `publication_XXXX` (unclear how different); This also could be the same as `development_date` above.
* `cellannotation_version` :fire::fire::fire: : **CAP** taxonomy annotation version; required by CAP; also publication_XXXX (unclear how different). I'm also not sure how this differs from the `cellannotation_schema_version`.
* `dataset_url`: file location, as defined above.
* `matrix_file_id`: file location, as defined above.
* `author_list`: list of taxonomy authors; see above
* `[additional information]` :fire::fire::fire: : ***Placeholder for several other (seemingly redundant) fields required by external tools (e.g., CAP, cellxgene) that I want to capture here. It may or may not make sense to spell them all out.***


END OF SCHEMA
