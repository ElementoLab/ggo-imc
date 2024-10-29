
import os
from glob import glob
from pathlib import Path
import tifffile
from tqdm import tqdm

import pandas as pd
import numpy as np

import scanpy as sc
import anndata

# from utag import utag
import seaborn as sns
import matplotlib.pyplot as plt
import imc_analysis as imc

import anndata
import warnings
warnings.simplefilter("ignore", UserWarning)

import matplotlib
sc.settings.set_figure_params(dpi=200, dpi_save=300, fontsize=12)
matplotlib.rcParams["pdf.fonttype"] = 42
matplotlib.rcParams["ps.fonttype"] = 42
matplotlib.rcParams["axes.grid"] = False
#matplotlib.use('Agg')


import yaml

metadata = load('metadata/ggo_config.yml')

adata_dict = dict()
utag_dict = dict()

for p in ['PANEL_G', 'PANEL_H']:
    if os.path.exists(metadata[p]['AnnData']['phenotyped_file_name']):

        adata_dict[p] = sc.read(metadata[p]['AnnData']['phenotyped_file_name'])
        utag_dict[p] = sc.read(metadata[p]['AnnData']['utag_file_name'])
    else:
        print(f"Reading {metadata[p]['AnnData']['phenotyped_file_name']}...")
        adata_dict[p] = sc.read(metadata[p]['AnnData']['phenotyped_file_name'])

        utag_results = utag(
            adata_dict[p],
            slide_key = "roi",
            max_dist = 20,
            normalization_mode = 'l1_norm',
            apply_clustering = False,
        )

        utag_results.write(metadata[p]['AnnData']['utag_file_name'])
        utag_dict[p] = utag_results

for p in ['PANEL_G', 'PANEL_H']:
    print(f'K means: {p}')
    from sklearn.cluster import DBSCAN, MiniBatchKMeans
    # clusters = DBSCAN(min_samples = 10000).fit(utag_dict[p].X)
    # # get cluster labels
    # utag_dict[p].obs['UTAG DB_scan'] = clusters.labels_
    mbk = MiniBatchKMeans(
        init ='k-means++',
        n_clusters = 20,
        random_state = 42,
    )
    mbk.fit(utag_dict[p].X)
    utag_dict[p].obs['UTAG K means'] = (mbk.labels_).astype(str)


cluster_id = 'UTAG K means'
for p in ['PANEL_G', 'PANEL_H']:

    metadata = load('metadata/ggo_config.yml')
    var_names = metadata[p]['var_celltype_groups']
    l = [set(l) for l in list(var_names.values())]
    varlist = list(set.union(*l))

    celltype_map = metadata[f'{p}_ue_20']
    celltype_map2 = {f'{k}: {celltype_map[k]}' for k in celltype_map}
    utag_dict[p].obs['uE'] = pd.Categorical(utag_dict[p].obs[cluster_id].replace(celltype_map))
    utag_dict[p].obs['mapped'] = pd.Categorical(utag_dict[p].obs[cluster_id].replace(celltype_map2))
    utag_dict[p].obs['uE_broad'] = pd.Categorical(utag_dict[p].obs['uE'].str.split(' (', regex=False).str[0])
    # utag_dict[p] = utag_dict[p][utag_dict[p].obs[c].isin(celltype_map.values()),:]


# adata_dict = dict()
# utag_dict = dict()

for p in ['PANEL_G', 'PANEL_H']:
    print(f"Reading {metadata[p]['AnnData']['phenotyped_umap_name']}...")
    # adata_dict[p] = sc.read(metadata[p]['AnnData']['phenotyped_umap_name'])
    adata_dict[p].obs['PID'] = adata_dict[p].obs['description'].str[:6]
    adata_dict[p].obs['radio'] = pd.Categorical(adata_dict[p].obs['radio'], categories = ['N', 'PNS', 'PS', 'S'], ordered=True)
    # adata_dict[p].obs['pathology'] = pd.Categorical(adata_dict[p].obs['pathology'], categories = ['Normal', 'AIS', 'MIA', 'IAC'], ordered = True)
    adata_dict[p].obs['pathology'] = adata_dict[p].obs['pathology'].astype(str).str.replace('Normal', 'N')
    adata_dict[p].obs['pathology'] = pd.Categorical(adata_dict[p].obs['pathology'], categories = ['N', 'AIS', 'MIA', 'IAC'], ordered = True)
    # utag_dict[p] = sc.read(metadata[p]['AnnData']['utag_labeled_name'])


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


import seaborn as sns
for p in ['PANEL_G', 'PANEL_H']:
    for celltype in ['celltype', 'celltype_broad']:
        for uE in ['uE', 'uE_broad']:
            count = adata_dict[p].obs.groupby([celltype, uE]).count()
            count = count / 1000
            g = sns.FacetGrid(
                data=count["sample"].reset_index(),
                hue=celltype,
                col=uE,
                col_wrap=2,
                aspect=1.3,
                palette="colorblind",
                sharex=False,
                height = 3,
            )

            g.map(sns.barplot, "sample",  celltype)

            # iterate over axes of FacetGrid
            for i, ax in enumerate(g.axes.flat):
                labels = ax.get_yticklabels()  # get x labels
                #ax.set_xticklabels(labels, rotation=90, fontsize=10)  # set new labels
                titles = ax.get_title()
                titles = titles.replace(f"{uE} = ", "")
                ax.set_title(titles, fontsize=12)
                ax.set_ylabel("", fontsize=12)
                ax.set_xlabel("")
            plt.tight_layout()
            plt.xlabel('x1000 cells')
            plt.savefig(f"figures/ue/{p}_{uE}_{celltype}_composition.pdf")

            # g = sns.FacetGrid(
            #     data=adata_dict[p].obs.groupby([celltype, "mapped"]).count()["sample"].reset_index(),
            #     hue=celltype,
            #     col="mapped",
            #     col_wrap=2,
            #     aspect=1.3,
            #     palette="colorblind",
            #     sharex=False,
            #     height = 6,
            # )

            # g.map(sns.barplot, "sample",  celltype)

            # # iterate over axes of FacetGrid
            # for i, ax in enumerate(g.axes.flat):
            #     labels = ax.get_yticklabels()  # get x labels
            #     #ax.set_xticklabels(labels, rotation=90, fontsize=10)  # set new labels
            #     titles = ax.get_title()
            #     titles = titles.replace(f"{uE} = ", "")
            #     ax.set_title(titles, fontsize=12)
            #     ax.set_ylabel("", fontsize=12)
            #     ax.set_xlabel("")
            # plt.tight_layout()
            # plt.savefig(f"figures/ue/{p}_{uE}_mapped_composition.pdf")

