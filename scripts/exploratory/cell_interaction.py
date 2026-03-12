import scanpy as sc
# import rapids_singlecell as rsc
import squidpy as sq
import anndata
from glob import glob
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import numpy as np
import yaml
import os
# import imc_analysis as imc
from pathlib import Path
sc.settings.set_figure_params(dpi=200, dpi_save=300, fontsize=12)

from tqdm import tqdm
import matplotlib
matplotlib.rcParams["pdf.fonttype"] = 42
matplotlib.rcParams["ps.fonttype"] = 42
matplotlib.rcParams["axes.grid"] = False
matplotlib.use('Agg')

from scripts.load_yaml import load
metadata = load('metadata/ggo_config.yml')

adata_dict = dict()

# read in data
for p in ['PANEL_G', 'PANEL_H']:
    adata_dict[p] = sc.read(metadata[p]['AnnData']['phenotyped_umap_name'])

# celltype_broad_map = {
#     'PANEL_G': {
#         'Epi.-like': 'Epithelial-like',
#         'Epi. Prol.': 'Epithelial-like (Ki67+)',
#         'Endo.': 'Endothelial',
#         'Fib.': 'Fibroblast',
#         'Mesen.-like': 'Mesenchymal-like',
#         'Mac.': 'Macrophage',
#         'NK': 'NK',
#         'T reg': 'T reg',
#         'CD4 T': 'CD4 T',
#         'CD8 T': 'CD8 T',
#         'B':'B',},
#     'PANEL_H': {
#         'Tumor-like': 'Tumor-like',
#         'Epi.-like': 'Epithelial-like',
#         'Endo.': 'Endothelial',
#         'Fib.': 'Fibroblast',
#         'Mac.': 'Macrophage',
#         'Mast': 'Mast',
#         'NK': 'NK',
#         'Neut.': 'Neutrophil',
#         'PMN-MDSC': 'PMN-MDSC',
#         'B': 'B',
#         'CD4 T': 'CD4 T',
#         'CD8 T':'CD8 T',}
# }

for p in ['PANEL_G', 'PANEL_H']:
    # adata_dict[p].obs['celltype_broad'] = adata_dict[p].obs['celltype_broad'].astype(str).replace(celltype_broad_map[p])
    # adata_dict[p].obs['celltype_broad'] = pd.Categorical(adata_dict[p].obs['celltype_broad'], categories = celltype_broad_map[p].values(), ordered = True)
    adata_dict[p] = adata_dict[p][adata_dict[p].obs['celltype_broad'].isin(adata_dict[p].obs['celltype_broad'].cat.categories)]
    adata_dict[p].obs['radio'] = pd.Categorical(adata_dict[p].obs['radio'], categories = ['N', 'PNS', 'PS', 'S'], ordered=True)
    adata_dict[p].obs['pathology'] = pd.Categorical(adata_dict[p].obs['pathology'], categories = ['Normal', 'AIS', 'MIA', 'IAC'], ordered = True)
 

# for p in ['PANEL_G', 'PANEL_H']:
#     adata_dict[p].write(metadata[p]['AnnData']['phenotyped_umap_name'])
    
interactions = dict()
for feature in ['radio', 'pathology']:
    interactions[feature] = {
        'PANEL_G': dict(),
        'PANEL_H': dict()
    }

# apply spatial neighborhood test
for p in ['PANEL_G', 'PANEL_H']:
    for feature in ['radio', 'pathology']:
        for f in tqdm(adata_dict[p].obs[feature].cat.categories):
        
            adata = adata_dict[p][adata_dict[p].obs[feature] == f].copy()
            # build graph
            sq.gr.spatial_neighbors(adata, coord_type="generic", spatial_key="spatial", library_key = 'roi', radius = 40)
            # test enrichment
            sq.gr.nhood_enrichment(adata, cluster_key="celltype_broad", show_progress_bar = False, seed = 0)
            # store as data frame of celltype by celltype
            ct = adata_dict[p].obs['celltype_broad'].cat.categories
            interactions[feature][p][f] = pd.DataFrame(adata.uns['celltype_broad_nhood_enrichment']['zscore'], index = ct, columns = ct)


# plot interaction heatmap
for feature in ['radio', 'pathology']:
    for p in ['PANEL_G', 'PANEL_H']:
        for f in interactions[feature][p]:
            sns.heatmap(interactions[feature][p][f], cmap = 'vlag', center = 0, vmax = 50, vmin = -50)
            plt.title(f'Interaction {feature}: {f}')
            os.makedirs(f'figures/{p}/interactions/{feature}/{f}/', exist_ok = True)
            plt.savefig(f'figures/{p}/interactions/{feature}/{f}/heatmap.pdf', bbox_inches = 'tight')
            plt.savefig(f'figures/{p}/interactions/{feature}/{f}/heatmap.png', bbox_inches = 'tight')
            plt.close()


# plot differential interaction heatmap and clustermap
feature = 'radio'
for p in ['PANEL_G', 'PANEL_H']:
    for f in interactions[feature][p]:
        os.makedirs(f'figures/{p}/differential_interaction/{feature}/{f}/', exist_ok = True)
        sns.heatmap(interactions[feature][p][f] - interactions[feature][p]['N'], cmap = 'vlag', center = 0, vmax = 50, vmin = -50)
        plt.title(f'Differential Interaction {feature}: {f}')
        plt.savefig(f'figures/{p}/differential_interaction/{feature}/{f}/heatmap.pdf', bbox_inches = 'tight')
        plt.savefig(f'figures/{p}/differential_interaction/{feature}/{f}/heatmap.png', bbox_inches = 'tight')
        plt.close()

        sns.clustermap(interactions[feature][p][f] - interactions[feature][p]['N'], cmap = 'vlag', center = 0, vmax = 50, vmin = -50)
        plt.title(f'Differential Interaction {feature}: {f}')
        plt.savefig(f'figures/{p}/differential_interaction/{feature}/{f}/clustmap.pdf', bbox_inches = 'tight')
        plt.savefig(f'figures/{p}/differential_interaction/{feature}/{f}/clustmap.png', bbox_inches = 'tight')
        plt.close()
    
feature = 'pathology'
for p in ['PANEL_G', 'PANEL_H']:
    for f in interactions[feature][p]:
        os.makedirs(f'figures/{p}/differential_interaction/{feature}/{f}/', exist_ok = True)
        sns.heatmap(interactions[feature][p][f] - interactions[feature][p]['Normal'], cmap = 'vlag', center = 0, vmax = 50, vmin = -50)
        plt.title(f'Differential Interaction {feature}: {f}')
        plt.savefig(f'figures/{p}/differential_interaction/{feature}/{f}/heatmap.pdf', bbox_inches = 'tight')
        plt.savefig(f'figures/{p}/differential_interaction/{feature}/{f}/heatmap.png', bbox_inches = 'tight')
        plt.close()

        sns.clustermap(interactions[feature][p][f] - interactions[feature][p]['Normal'], cmap = 'vlag', center = 0, vmax = 50, vmin = -50)
        plt.title(f'Differential Interaction {feature}: {f}')
        plt.savefig(f'figures/{p}/differential_interaction/{feature}/{f}/clustmap.pdf', bbox_inches = 'tight')
        plt.savefig(f'figures/{p}/differential_interaction/{feature}/{f}/clustmap.png', bbox_inches = 'tight')
        plt.close()
