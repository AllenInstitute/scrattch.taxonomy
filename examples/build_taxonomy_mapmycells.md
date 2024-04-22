# Tutorial: Building a MapMyCells taxonomy 

This tutorial shows how to set up a taxonomy using the python version of cell type mapper (MapMyCells) for running HANN python mapping algorithm against.

## Overview
#### Required inputs:

* Count matrix (gene x cell), with genes as colnames and sample identifiers as rownames

#### Additional prerequisites:

* Reference dataset, this tutorial uses Siletti subsampled dataset:
  /allen/programs/celltypes/workgroups/hct/cellTaxonomy/adult-human-brain_v1/additional_files/processed/CSR/human_whole_brain_subcluster_centroid_subsampled.h5ad
* cell_type_mapper library https://github.com/AllenInstitute/cell_type_mapper


## How to build the taxonomy
### Step1: Python environment installation

*Note*: skip the substeps 1 and 2 if you already have a python environment set up. 

1. Install miniconda:

```
mkdir -p ~/miniconda3
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda3/miniconda.sh
bash ~/miniconda3/miniconda.sh -b -u -p ~/miniconda3
rm -rf ~/miniconda3/miniconda.sh
```

*Note*:
* These commands will install the miniconda to your home directory. To install it in a different directory, replace '~' with the path to the folder in all of the command below. 
* For more help with installing miniconda visit [https://docs.anaconda.com/free/miniconda/miniconda-install.html](https://docs.anaconda.com/free/miniconda/)

2. Create a python environment: 
```
conda create --name <your-env-name> -y
```

More information: https://conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html

*Note*:
The conda enviromnet needs to be activated before running cell_type_mapper

3. Install cell_type_mapper library:

```
conda activate <your-env-name>
git clone https://github.com/AllenInstitute/cell_type_mapper.git <path_to_your_directory>
cd <path_to_your_directory>/cell_type_mapper
pip install -r requirements.txt
pip install -e .
```

*The coding steps above are*:
* Activate your conda environment.
* Clone cell_type_mapper repository to your local directory.
* Go to cell_type_mapper folder in your local directory.
* Install the required packages (dependencies).
* Install the cell_type_mapper package itself. Note: this needs to be run from the root directory of this repository, i.e. <path_to_your_directory>/cell_type_mapper

### Step 2: Build the taxonomy

*Note*: conda environment needs to be activated to run the code below. Run `conda activate <your-env-name>` to activate.
  
1. Build precompute_stats file:

*<path_to_your_taxonomy_files_directory> below is a full path to an already created folder where the build files are saved.* 

*temp folder is also needs to be created and passed in --tmp_dir*

```
python -m cell_type_mapper.cli.precompute_stats_scrattch \
--h5ad_path /allen/programs/celltypes/workgroups/hct/cellTaxonomy/adult-human-brain_v1/additional_files/processed/CSR/human_whole_brain_subcluster_centroid_subsampled.h5ad \
--hierarchy '["supercluster_term", "cluster_id", "subcluster_id"]' \
--output_path <path_to_your_taxonomy_files_directory>/siletti_hmba_subsampled_precompute_stats.h5 \
--normalization raw \
--tmp_dir <path_to_your_taxonomy_files_directory>/temp/
```


*Note*:
* Create a folder in your own directory where the build taxonomy files are saved to later be used for mapping and pass instead of <path_to_your_taxonomy_files_directory>.
* Change the query h5ad_path parameter above with any other h5ad taxonomy file location.
* Change the name of the precompute_stats file (siletti_hmba_subsampled_precompute_stats.h5) if you're using a different taxonomy.
* Change the hierarchy based on the hierarchy of the taxonomy; for example '["class", "subclass", "cluster"]'.
* Change normalization to either 'raw' or 'log2CPM' based on the count matrix of the taxonomy.
* Create a temp folder in either the directory where build taxonomy files are saved, or any other directory. This folder is just for the temporary files that the algorithm creates and later deletes.
  
  
2. Build reference_markers file:
```
python -m cell_type_mapper.cli.reference_markers \
--precomputed_path_list "['<path_to_your_taxonomy_files_directory>/siletti_hmba_subsampled_precompute_stats.h5']" \
--output_dir <path_to_your_taxonomy_files_directory> \
--tmp_dir <path_to_your_taxonomy_files_directory>/temp/
```


*Note*:
* Pass a different --precomputed_path_list parameter if it was changed in the previous step.
* Use same <path_to_your_taxonomy_files_directory> and temp folder directory as in the previous step.

  
3. Build query_markers file
```
python -m cell_type_mapper.cli.query_markers \
--reference_marker_path_list '["<path_to_your_taxonomy_files_directory>/reference_markers.h5"]' \
--output_path <path_to_your_taxonomy_files_directory>/siletti_hmba_subsampled_query_markers.json
```


*Note*:
* Change the name of the query_markers file (siletti_hmba_subsampled_query_markers.json) if you're using a different taxonomy.
* Use same <path_to_your_taxonomy_files_directory> as in the previous steps.
  
