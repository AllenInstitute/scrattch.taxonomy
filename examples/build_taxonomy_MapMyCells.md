# Tutorial: Building a MapMyCells taxonomy

This tutorial shows how to set up a taxonomy using the python version of cell type mapper (MapMyCells) for running HANN python mapping algorithm against.

## Overview
#### Required inputs:

* Count matrix (cell x gene), with genes as colnames and sample identifiers as rownames.

#### Additional prerequisites:

* Reference taxonomy in .h5ad format with a predefined heirarchy of cell annotations.
* cell_type_mapper library https://github.com/AllenInstitute/cell_type_mapper

## How to build the taxonomy

### Step1: Python environment installation [Link](https://github.com/AllenInstitute/scrattch.mapping/blob/inkar-HANN-python-tutorial/examples/setup_env_MapMyCells.md)

### Step 2: Build the taxonomy

*Note*: conda environment needs to be activated to run the code below. Run `conda activate <your-env-name>` to activate.
  
#### 2.1. Build precompute_stats file:

*<path_to_your_taxonomy_files_directory> below is a full path to an already created folder where the build files are going to be saved.* 

*temp folder is also needs to be created and passed as --tmp_dir parameter*

```
python -m cell_type_mapper.cli.precompute_stats_scrattch \
--h5ad_path <path_to_your_taxonomy_files>/<your_taxonomy>.h5ad \
--hierarchy '["class", "subclass", "cluster"]' \
--output_path <path_to_your_taxonomy_files_directory>/<your_taxonomy>_precompute_stats.h5 \
--normalization raw \
--tmp_dir <path_to_your_taxonomy_files_directory>/temp/
```

*Note*:
* Create a folder in your own directory where the built taxonomy files are saved, to later be used for mapping, and pass instead of <path_to_your_taxonomy_files_directory>.
* Change the reference taxonomy file --h5ad_path parameter above to any other h5ad taxonomy file location.
* Change the name of the precompute_stats file based on your taxonomy.
* Change the hierarchy based on the hierarchy of your taxonomy; for example '["class", "subclass", "cluster"]'.
* Change normalization to either 'log2CPM' or keep as 'raw' based on the count matrix of your taxonomy.
* Create a temp folder in either the directory where build taxonomy files are saved, or any other directory. This folder is just for the temporary files that the algorithm creates and later deletes.

For more information about other command line parameters, run in your terminal:
```
python -m cell_type_mapper.cli.precompute_stats_scrattch --help
```
  
#### 2.2. Build reference_markers file:
```
python -m cell_type_mapper.cli.reference_markers \
--precomputed_path_list "['<path_to_your_taxonomy_files_directory>/<your_taxonomy>_precompute_stats.h5']" \
--output_dir <path_to_your_taxonomy_files_directory> \
--tmp_dir <path_to_your_taxonomy_files_directory>/temp/
```

*Note*:
* Pass a different --precomputed_path_list parameter if it was changed in the previous step.
* Use the same <path_to_your_taxonomy_files_directory> and temp folder directory as in the previous step.

For more information about other command line parameters, run in your terminal:
```
python -m cell_type_mapper.cli.reference_markers --help
```
  
#### 2.3. Build query_markers file
```
python -m cell_type_mapper.cli.query_markers \
--reference_marker_path_list '["<path_to_your_taxonomy_files_directory>/reference_markers.h5"]' \
--output_path <path_to_your_taxonomy_files_directory>/<your_taxonomy>_query_markers.json
```

*Note*:
* Change the name of the query_markers file if you're using a different taxonomy.
* Use the same <path_to_your_taxonomy_files_directory> as in the previous steps.
* The reference_markers file is not needed for mapping once the query_markers file is generated. 

For more information about other command line parameters, run in your terminal:
```
python -m cell_type_mapper.cli.query_markers --help
```
