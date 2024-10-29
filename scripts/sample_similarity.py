
import os
from glob import glob
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
from scripts.load_yaml import load

import matplotlib
sc.settings.set_figure_params(dpi=200, dpi_save=300, fontsize=12)
matplotlib.rcParams["pdf.fonttype"] = 42
matplotlib.rcParams["ps.fonttype"] = 42
matplotlib.rcParams["axes.grid"] = False
matplotlib.use('Agg')

# load data
metadata = load('metadata/ggo_config.yml')

# Read in data
adata_dict = dict()
for p in ['PANEL_G', 'PANEL_H']:
    print(f"Reading {metadata[p]['AnnData']['phenotyped_umap_name']}...")
    adata_dict[p] = sc.read(metadata[p]['AnnData']['phenotyped_umap_name'])
    adata_dict[p].obs['PID'] = adata_dict[p].obs['description'].str[:6]
    adata_dict[p].obs['radio'] = pd.Categorical(adata_dict[p].obs['radio'], categories = ['N', 'PNS', 'PS', 'S'], ordered=True)
    adata_dict[p].obs['pathology'] = pd.Categorical(adata_dict[p].obs['pathology'], categories = ['Normal', 'AIS', 'MIA', 'IAC'], ordered = True)
    del adata_dict[p].obs['ROI_area']

conditions = metadata['CONDITIONS']

# Get cell type density matrix
celltype = 'celltype_broad'
density_dict = dict()
for p in ['PANEL_G', 'PANEL_H']:
    adata = adata_dict[p]
    adata.obs['all_cells'] = 'all_cells'
    
    density_dict[p] = imc.tl.celltype_density(
        adata,
        celltype = celltype,
        condition_keys = conditions + ['PID', 'GGO ID'])
    if 'GGO ID' in density_dict[p].obs:
        idx_keep = density_dict[p].obs[['GGO ID']].drop_duplicates().index
        density_dict[p] = density_dict[p][idx_keep]

        density_dict[p].obs.set_index('GGO ID', inplace = True)

# Get common rois and celltypes
common_roi = pd.Index.intersection(*[density_dict[p].obs.index for p in density_dict])
common_celltype = pd.Index.intersection(*[density_dict[p].var.index for p in density_dict])
common_celltype = common_celltype.drop('Low Expr.')

# Save panels
for p in density_dict:
    density_dict[p] = density_dict[p][common_roi,common_celltype]
    density_dict[p].obs['panel'] = p

# concatenate densities
density = anndata.concat([density_dict[p] for p in density_dict])
density.obs['sample'] = pd.Categorical(density.obs['GGO ID'].str[6:], categories = ['N', 'T_1', 'T_2'])

# plot similarity per celltype
data = density.to_df()
pal = sns.color_palette(['lightgray'], density.shape[1])

# swarm box line plot
fig, axes = plt.subplots(2,3,figsize = (6,4), dpi = 300)
for i, ax in enumerate(axes.flatten()):
    if i < len(density.var.index):
        celltype = density.var.index[i]
        sns.boxplot(x="sample", y=data[celltype], data=density.obs, fill = False, fliersize = 0, color ='.25', ax = ax)
        sns.swarmplot(x="sample", y=data[celltype], data=density.obs, color=".25", size = 2, ax = ax)
        sns.lineplot(x="sample", y=data[celltype], hue='PID', data=density.obs,
                     estimator=None, legend=False, palette = pal, alpha = 0.7, size = 0.5, ax = ax)
        ax.set_xlabel('')
        ax.set_ylabel('')
        ax.set_title(celltype)
    else:
        ax.remove()
plt.tight_layout()
plt.savefig('figures/sample_similarity_boxswarmlineplot.pdf')
plt.close()

# plot similarity per panel per celltype as scatter
import scipy as sp
fig, axes = plt.subplots(2,3,figsize = (6,4), dpi = 300)
for i, ax in enumerate(axes.flatten()):
    if i < len(density.var.index):
        celltype = density.var.index[i]
        d_ = density.obs.copy()
        d_[celltype] = data[celltype]
        d_ = d_.pivot_table(index = 'PID', columns = 'sample', values = celltype)
        arr = d_[['T_1','T_2']].dropna()
        arr['T1'] = arr[['T_1', 'T_2']].min(axis = 1)
        arr['T2'] = arr[['T_1', 'T_2']].max(axis = 1)
        
        # arr = np.log10(1+arr)
        sns.regplot(x="T1", y="T2", data=arr, ax = ax, color=".3", truncate = False, line_kws=dict(color="r", linewidth = 0.5, alpha = 0.5), scatter_kws={'s':0.5})
        ax.set_title(celltype)
        ax.set_xlabel('')
        ax.set_ylabel('')
        # ax.set_xscale('log')
        # ax.set_yscale('log')
        # ax.set_xlim(1.5,3.5)
        # ax.set_ylim(1.5,3.5)
        r, p_ = sp.stats.spearmanr(arr['T1'], arr['T2'])
        ax.text(.05, .8, 'r={:.2f}, p={:.2g}'.format(r, p_),
                transform=ax.transAxes)

    else:
        ax.remove()
sns.despine()
plt.tight_layout()
plt.savefig('figures/sample_similarity_regplot_log.pdf')
plt.close()

# save similarity per patient as heatmap input matrix
heatmap_data = density.copy()
sample_groups = ['N', 'T_1', 'T_2']
panels = ['PANEL_G', 'PANEL_H']

fractions = []
for sample in sample_groups:
    for panel in panels:
        # for celltype in heatmap_data.var.index:
        fraction = heatmap_data[
            (heatmap_data.obs['sample'] == sample) &
            (heatmap_data.obs['panel'] == panel)
        ].copy()
        fraction.var['sample'] = sample
        fraction.var['panel'] = panel
        fractions.append(fraction)

common_pid = list(set.intersection(*[set(f.obs['PID']) for f in fractions]))
fractions2 = []
for f in fractions:
    f2 = f[f.obs['PID'].isin(common_pid)].copy()
    f2.obs.set_index('PID',inplace=True)
    fractions2.append(f2)

heatmap_adata = anndata.concat(fractions2, axis = 1)
heatmap_adata.var.reset_index(inplace = True)
heatmap_adata.obs.reset_index(inplace = True)

p = sns.clustermap(heatmap_adata.to_df(), cmap='RdBu_r', linewidths=0.5)
order = p.dendrogram_row.reordered_ind

heatmap_adata = heatmap_adata[heatmap_adata.obs.loc[order].index, heatmap_adata.var.sort_values(['celltype_broad', 'sample', 'panel']).index]
heatmap_adata.write('results/cell heatmap density.h5ad')