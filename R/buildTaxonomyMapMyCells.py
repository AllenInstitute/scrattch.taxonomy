# anndata_path="/allen/programs/celltypes/workgroups/rnaseqanalysis/EvoGen/Team/Inkar/taxonomies/tasic/Tasic2016.h5ad"
# hierarchy = ["broad_type_label", "primary_type_label"]
#hierarchy = ["supercluster_term", "cluster_id", "subcluster_id"]
# n_processors = 3
# normalization = 'log2CPM'
from cell_type_mapper.cli.precompute_stats_scrattch import (PrecomputationScrattchRunner)
from cell_type_mapper.cli.reference_markers import (ReferenceMarkerRunner)
from cell_type_mapper.cli.query_markers import (QueryMarkerRunner)
from cell_type_mapper.utils import output_utils 
from cell_type_mapper.utils import anndata_utils

import time
import os
import json

def buildMapMyCellsTaxonomy(anndata_path, hierarchy, tmp_dir, n_processors, normalization):
    # add try catch, save AIT.anndata? or not, redo the imports
    delete_tmp_folder = False
    if not tmp_dir:
        tmp_dir = "./temp_folder_" + time.strftime("%Y%m%d-%H%M%S")
        delete_tmp_folder = True

    precomp_stats_filename = "precompute_stats_" + time.strftime("%Y%m%d-%H%M%S") + ".h5"
    precomp_stats_output_path = os.path.join(tmp_dir, precomp_stats_filename)

    # PRECOMPUTE STATS
    precomp_stats_config = {
        'h5ad_path': anndata_path,
        'n_processors': n_processors,
        'normalization': normalization,
        'tmp_dir': tmp_dir,
        'output_path': precomp_stats_output_path,
        'hierarchy': hierarchy
    }

    precomp_stats_runner = PrecomputationScrattchRunner(
        args=[],
        input_data=precomp_stats_config)
    precomp_stats_runner.run()

    # ADD TO UNS
    output_utils.precomputed_stats_to_uns(
        precomputed_stats_path=precomp_stats_output_path, 
        h5ad_path=anndata_path, 
        uns_key="MapMyCells_precomp_stats"
    )

    # REFERENCE MARKERS
    ref_markers_config = {
        'n_processors': n_processors,
        'precomputed_path_list': [precomp_stats_output_path],
        'output_dir': tmp_dir,
        'tmp_dir': tmp_dir
    }

    ref_markers_runner = ReferenceMarkerRunner(
        args=[], 
        input_data=ref_markers_config)
    ref_markers_runner.run()

    # QUERY MARKERS
    ref_markers_file_path = os.path.join(tmp_dir, "reference_markers.h5")
    query_markers_filename = "query_markers_" + time.strftime("%Y%m%d-%H%M%S") + ".h5"
    query_markers_output_path = os.path.join(tmp_dir, query_markers_filename)

    query_markers_config = {
        'query_path': anndata_path,
        'reference_marker_path_list': [ref_markers_file_path],
        'n_processors': n_processors,
        'output_path': query_markers_output_path,
        'tmp_dir': tmp_dir}

    query_markers_runner = QueryMarkerRunner(
        args=[],
        input_data=query_markers_config)
    query_markers_runner.run()

    query_markers_data = json.load(open(query_markers_output_path))
    query_markers_uns_key = "MapMyCells_query_markers"
    anndata_utils.update_uns(
        h5ad_path=anndata_path,
        new_uns={query_markers_uns_key: query_markers_data},
        clobber=False
    )

    # redo for if execution haults before reaching this line
    os.remove(precomp_stats_output_path)
    os.remove(ref_markers_file_path)
    os.remove(query_markers_output_path)
    if delete_tmp_folder:
        os.rmdir(tmp_dir)

    print("Finished building MapMyCells taxonomy.")