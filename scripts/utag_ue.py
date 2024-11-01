
import pandas as pd
import numpy as np

import scanpy as sc
import anndata

# from utag import utag
import seaborn as sns
import matplotlib.pyplot as plt
import imc_analysis as imc

metadata = imc.utils.parse_yaml('metadata/ggo_config.yml')


# # Method for UTAG #
# import os

# adata_dict = dict()
# utag_dict = dict()
# for p in ['PANEL_G', 'PANEL_H']:

#     if os.path.exists(metadata[p]['AnnData']['phenotyped_file_name']):

#         adata_dict[p] = sc.read(metadata[p]['AnnData']['phenotyped_file_name'])
#         utag_dict[p] = sc.read(metadata[p]['AnnData']['utag_file_name'])
#     else:
#         print(f"Reading {metadata[p]['AnnData']['phenotyped_file_name']}...")
#         adata_dict[p] = sc.read(metadata[p]['AnnData']['phenotyped_file_name'])

#         utag_results = utag(
#             adata_dict[p],
#             slide_key = "roi",
#             max_dist = 20,
#             normalization_mode = 'l1_norm',
#             apply_clustering = False,
#         )

#         utag_results.write(metadata[p]['AnnData']['utag_file_name'])
#         utag_dict[p] = utag_results

# for p in ['PANEL_G', 'PANEL_H']:
#     print(f'K means: {p}')
#     from sklearn.cluster import DBSCAN, MiniBatchKMeans
#     # clusters = DBSCAN(min_samples = 10000).fit(utag_dict[p].X)
#     # # get cluster labels
#     # utag_dict[p].obs['UTAG DB_scan'] = clusters.labels_
#     mbk = MiniBatchKMeans(
#         init ='k-means++',
#         n_clusters = 20,
#         random_state = 42,
#     )
#     mbk.fit(utag_dict[p].X)
#     utag_dict[p].obs['UTAG K means'] = (mbk.labels_).astype(str)


# cluster_id = 'UTAG K means'
# for p in ['PANEL_G', 'PANEL_H']:

#     metadata = imc.utils.parse_yaml('metadata/ggo_config.yml')
#     var_names = metadata[p]['var_celltype_groups']
#     l = [set(l) for l in list(var_names.values())]
#     varlist = list(set.union(*l))

#     celltype_map = metadata[f'{p}_ue_20']
#     celltype_map2 = {f'{k}: {celltype_map[k]}' for k in celltype_map}
#     utag_dict[p].obs['uE'] = pd.Categorical(utag_dict[p].obs[cluster_id].replace(celltype_map))
#     utag_dict[p].obs['mapped'] = pd.Categorical(utag_dict[p].obs[cluster_id].replace(celltype_map2))
#     utag_dict[p].obs['uE_broad'] = pd.Categorical(utag_dict[p].obs['uE'].str.split(' (', regex=False).str[0])
#     # utag_dict[p] = utag_dict[p][utag_dict[p].obs[c].isin(celltype_map.values()),:]
#     # utag_dict[p].write(metadata[p]['AnnData']['utag_labeled_name'])

adata_dict = dict()
utag_dict = dict()
for p in ['PANEL_G', 'PANEL_H']:
    print(f"Reading {metadata[p]['AnnData']['phenotyped_umap_name']}...")
    adata_dict[p] = sc.read(metadata[p]['AnnData']['phenotyped_umap_name'])
    utag_dict[p] = sc.read(metadata[p]['AnnData']['utag_labeled_name'])

for p in adata_dict:
    if 'uE' in adata_dict[p].obs:
        del adata_dict[p].obs['uE']

    if 'uE_broad' in adata_dict[p].obs:
        del adata_dict[p].obs['uE_broad']

    obs = adata_dict[p].obs
    obs = obs.merge(utag_dict[p].obs[['uE','uE_broad']], left_index = True, right_index = True, how = 'left')
    adata_dict[p].obs['uE'] = obs['uE']
    adata_dict[p].obs['uE_broad'] = obs['uE_broad']
    adata_dict[p].write(metadata[p]['AnnData']['utag_labeled_name'])
