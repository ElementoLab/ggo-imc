
import os
from pathlib import Path
import tifffile
from tqdm import tqdm

import pandas as pd
import numpy as np

import scanpy as sc
import anndata

import seaborn as sns
import matplotlib.pyplot as plt

import anndata
import warnings
warnings.simplefilter("ignore", UserWarning)

import imc_analysis as imc

import matplotlib
sc.settings.set_figure_params(dpi=200, dpi_save=300, fontsize=20)
matplotlib.rcParams["pdf.fonttype"] = 42
matplotlib.rcParams["ps.fonttype"] = 42
matplotlib.rcParams["axes.grid"] = False
matplotlib.use('Agg')

from scripts.load_yaml import load
metadata = load('metadata/ggo_config.yml')


adata_dict = dict()
for p in ['PANEL_G', 'PANEL_H']:
    print(f"Reading {metadata[p]['AnnData']['phenotyped_umap_name']}...")
    adata_dict[p] = sc.read(metadata[p]['AnnData']['phenotyped_umap_name'])
    adata_dict[p].obs['PID'] = adata_dict[p].obs['description'].str[:6]
    adata_dict[p].obs['radio'] = pd.Categorical(adata_dict[p].obs['radio'], categories = ['N', 'PNS', 'PS', 'S'], ordered=True)
    adata_dict[p].obs['pathology'] = pd.Categorical(adata_dict[p].obs['pathology'], categories = ['Normal', 'AIS', 'MIA', 'IAC'], ordered = True)

celltypes = ['celltype', 'celltype_broad']

# for celltype in celltypes:
#     adata_dict[p] = adata_dict[p][adata_dict[p].obs[celltype] != 'Low Expr.']

print('Plot celltype heatmap')
for p in adata_dict:
    imc.pl.celltype_heatmap(
        adata_dict[p],
        cluster_ids = celltypes,
        var_names = metadata[p]['var_celltype_groups'],
        out_dir = f'figures/{p}/celltype/'
    )


print('Plot UMAP')
for p in adata_dict:
    adata = sc.pp.subsample(
        adata_dict[p],
        copy = True,
        n_obs = 100000,
        random_state = 0,
    )
    for celltype in metadata['CELLTYPES']:
        fig, ax = plt.subplots(1,1,figsize = (6,5), dpi = 300)
        sc.pl.umap(adata, color = celltype, frameon = False, show = False)

        os.makedirs(f'figures/{p}/{celltype}/', exist_ok = True)
        plt.tight_layout()
        plt.savefig(f'figures/{p}/{celltype}/umap.pdf')
        plt.close()
