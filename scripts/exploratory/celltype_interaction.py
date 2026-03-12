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
import networkx as nx
import warnings
warnings.simplefilter("ignore", UserWarning)

import matplotlib
matplotlib.rcParams["pdf.fonttype"] = 42
matplotlib.rcParams["ps.fonttype"] = 42
matplotlib.rcParams["axes.grid"] = False
matplotlib.use('Agg')

import yaml
metadata_filename = 'metadata/ggo_config.yml'

with open(metadata_filename, "r") as stream:
    try:
        metadata = yaml.safe_load(stream)
    except yaml.YAMLError as exc:
        print(exc)

adata_dict = dict()

print(f"Reading {metadata['PANEL_G']['AnnData']['phenotyped_file_name']}...")
adata_dict['PANEL_G'] = sc.read(metadata['PANEL_G']['AnnData']['phenotyped_file_name'])

print(f"Reading {metadata['PANEL_H']['AnnData']['phenotyped_file_name']}...")
adata_dict['PANEL_H'] = sc.read(metadata['PANEL_H']['AnnData']['phenotyped_file_name'])


def cell_type_colocalization(
    adata: anndata.AnnData,
    celltype_key: str = "celltype",
    max_dist: int = 40,
    n_iterations: int = 1000,
):
    import squidpy as sq
    a_ = adata.copy()
    sq.gr.spatial_neighbors(a_, radius=max_dist, coord_type="generic")

    G = nx.from_scipy_sparse_array(a_.obsp["spatial_connectivities"])

    utag_map = {i: x for i, x in enumerate(adata.obs[celltype_key])}
    nx.set_node_attributes(G, utag_map, name=celltype_key)

    adj, order = nx.linalg.attrmatrix.attr_matrix(G, node_attr=celltype_key)
    order = pd.Series(order).astype(adata.obs[celltype_key].dtype)
    freqs = pd.DataFrame(adj, order, order).fillna(0) + 1

    norm_freqs = correct_interaction_background_random(G, freqs, celltype_key, n_iterations)
    return norm_freqs


def correct_interaction_background_random(
    graph: nx.Graph,
    freqs: pd.DataFrame,
    attribute: str,
    n_iterations: int = 1000
):
    values = {x: graph.nodes[x][attribute] for x in graph.nodes}
    shuffled_freqs = list()
    for _ in range(n_iterations):
        g2 = graph.copy()
        shuffled_attr = pd.Series(values).sample(frac=1)
        shuffled_attr.index = values
        nx.set_node_attributes(g2, shuffled_attr.to_dict(), name=attribute)
        rf, rl = nx.linalg.attrmatrix.attr_matrix(g2, node_attr=attribute)
        rl = pd.Series(rl, dtype=freqs.index.dtype)
        shuffled_freqs.append(pd.DataFrame(rf, index=rl, columns=rl))#.fillna(0) + 1)
    shuffled_freq = pd.concat(shuffled_freqs)
    shuffled_freq = shuffled_freq.groupby(level=0).sum()
    shuffled_freq = shuffled_freq.fillna(0)
    
    fl = np.log(1+(freqs / freqs.values.sum()))
    sl = np.log(1+(shuffled_freq / shuffled_freq.values.sum()))
    # make sure both contain all edges/nodes
    fl = fl.reindex(sl.index, axis=0).reindex(sl.index, axis=1)
    sl = sl.reindex(fl.index, axis=0).reindex(fl.index, axis=1)
    return fl - sl

import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)

for panel in adata_dict:
    print(f'Calculating {panel} celltype interactions ...')
    adata = adata_dict[panel]
    for celltype in metadata['CELLTYPES']:
        for roi in tqdm(adata.obs['roi'].unique()):
            a = adata[adata.obs['roi'] == roi].copy()
            path = f'metadata/interactions/{panel}_{celltype}_{roi}.csv'
            if not os.path.exists(path):
                freq = cell_type_colocalization(a, celltype_key = celltype, n_iterations = 100)
                freq.to_csv(path)

celltype = "celltype_broad"
slide_key = "roi"
cond = 'radio'
from collections import defaultdict

from collections import defaultdict
for panel in adata_dict:
    for celltype in metadata['CELLTYPES']:
        for cond in metadata['CONDITIONS']:
            
            print(f'Calculating {panel} {celltype} interactions: {cond}')

            interactions = defaultdict(lambda: [])
            for c_ in tqdm(adata_dict[panel].obs[cond].cat.categories):
                
                adata = adata_dict[panel][adata_dict[panel].obs[cond] == c_].copy()
                for roi in adata.obs[slide_key].unique():
                    a = adata[adata.obs[slide_key] == roi].copy()
                    path = f'metadata/interactions/{panel}_{celltype}_{roi}.csv'
                    if os.path.exists(path):
                        freq = pd.read_csv(path, index_col = 0)
                    else:
                        freq = cell_type_colocalization(a, celltype_key = celltype, n_iterations = 100)
                        freq.to_csv(path)
                    df = freq.melt(ignore_index = False, value_name = 'log likelihood').reset_index()
                    df[slide_key] = roi
                    interactions[c_].append(df)

                interaction = pd.concat(interactions[c_])
                interaction.to_csv(f'metadata/interactions/{panel}_{celltype}_{cond}_{c_}.csv')
                mean_interaction = interaction.groupby(['index','variable']).mean().reset_index().pivot(index = 'index', columns = 'variable')
                mean_interaction.columns = mean_interaction.columns.get_level_values(1)
                mean_interaction.to_csv(f'metadata/interactions/mean_interaction_{panel}_{celltype}_{cond}_{c_}.csv')

                sns.heatmap(mean_interaction, cmap = 'vlag', center = 0)#, vmax = 1, vmin = -1)
                plt.title(f'Interaction {celltype} {cond}: {c_}')
                plt.savefig(f'figures/interactions/interaction_{panel}_{celltype}_{cond}_{c_}.pdf', bbox_inches = 'tight')
                plt.savefig(f'figures/interactions/interaction_{panel}_{celltype}_{cond}_{c_}.png', bbox_inches = 'tight')
                plt.close()


cond = 'pathology'
#cond = 'radio'

celltype = 'celltype_broad'
for celltype in metadata['CELLTYPES']:

    for panel in adata_dict:
        for cond in ['pathology', 'radio']:
            if cond == 'radio':
                N = 'N'
            elif cond == 'pathology':
                N = 'Normal'
            norm_interaction = pd.read_csv(f'metadata/interactions/mean_interaction_{panel}_{celltype}_{cond}_{N}.csv', index_col = 0)
            
            cat = adata_dict[panel].obs[cond].cat.categories
            for c_ in cat:
                interaction = pd.read_csv(f'metadata/interactions/mean_interaction_{panel}_{celltype}_{cond}_{c_}.csv', index_col = 0)
                sns.heatmap(
                    (interaction - norm_interaction).fillna(0),
                    xticklabels = True,
                    yticklabels = True,
                    cmap = 'vlag',
                    center = 0,
                    # vmax = 0.01,
                    # vmin = -0.01
                )
                plt.title(f'Interaction {panel.upper()} {celltype} {cond}: {c_} - Normal')
                plt.savefig(f'figures/interactions/{panel}_{celltype}_{cond}_{c_}-Normal_heatmap.pdf', bbox_inches = 'tight')
                plt.savefig(f'figures/interactions/{panel}_{celltype}_{cond}_{c_}-Normal_heatmap.png', bbox_inches = 'tight')
                plt.close()

                try:
                    sns.clustermap(
                        (interaction - norm_interaction).fillna(0),
                        xticklabels = True,
                        yticklabels = True,
                        cmap = 'vlag',
                        center = 0,
                        vmax = 0.005,
                        vmin = -0.005
                    )
                    plt.title(f'Interaction {panel.upper()} {celltype} {cond}: {c_} - Normal')
                    plt.savefig(f'figures/interactions/{panel}_{celltype}_{cond}_{c_}-Normal_clustermap.pdf', bbox_inches = 'tight')
                    plt.savefig(f'figures/interactions/{panel}_{celltype}_{cond}_{c_}-Normal_clustermap.png', bbox_inches = 'tight')
                    plt.close()
                except ValueError:
                    print('contains non-finite distance values')

            '''
            pattern = f'metadata/interactions/interaction_{panel}_{celltype}*_{cond}*'
            files = glob(pattern)
            tmp = []
            for file in files:
                df = pd.read_csv(file, index_col = 0)
                df[f'{cond}'] = file.split('_')[-1].replace('.csv', '')
                tmp.append(df)
            
            df = pd.concat(tmp)

            for i in ['index', 'variable', f'{cond}']:
                df[i] = pd.Categorical(df[i])
            l = len(df['index'].cat.categories)

            fig, axes = plt.subplots(l,l, dpi = 300, figsize = (l,l))

            for i, ax in tqdm(enumerate(axes.flatten())):
                i = int(i/l)
                v = i%l

                index = df['index'].cat.categories[i]
                variable = df['variable'].cat.categories[v]

                sns.boxplot(
                    data = df[(df['index'] == index) & (df['variable'] == variable)],
                    x = f'{cond}',
                    y = 'log likelihood',
                    palette = metadata[f'{cond}_color'],
                    ax = ax
                )
                ax.set_xlabel(variable)
                ax.set_ylabel(index)
            plt.suptitle(f'{cond} {celltype} Interaction Panel {panel}')
            plt.savefig(f'figures/{cond} {celltype} {panel}.pdf', bbox_inches = 'tight')
            plt.savefig(f'figures/{cond} {celltype} {panel}.png', bbox_inches = 'tight')
            plt.close()
            '''
