import scanpy as sc
import rapids_singlecell as rsc
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


def cluster(
    adata: anndata.AnnData,
    key_added: str = 'rep'
):
    a = adata.copy()

    # rsc.get.anndata_to_GPU(a)
    rsc.pp.pca(a)
    #rsc.pp.neighbors(a)
    sc.pp.neighbors(a, method='rapids')
    rsc.tl.leiden(a, key_added = key_added)
    # rsc.get.anndata_to_CPU(a)

    return a

def cluster_rep(
    adata: anndata.AnnData,
    n_reps: int = 3,
):
    
    adata = adata.copy()
    adata.obs['rep0'] = '1'
    a_list = [adata]
    a_list_next_rep = []
    for i in tqdm(range(n_reps)):
        
        for a_tmp in a_list:
            a_tmp2 = cluster(a_tmp, f'rep{i}')

            for c in a_tmp2.obs[f'rep{i}'].unique():
                a_list_next_rep.append(a_tmp2[a_tmp2.obs[f'rep{i}'] == c])

        a_list = a_list_next_rep
        a_list_next_rep = []

    a_clustered = anndata.concat(a_list)

    adata.obs[f'rep0'] = a_clustered.obs[f'rep0']
    for i in range(1, n_reps):
        adata.obs[f'rep{i}'] = adata.obs[f'rep{i-1}'].astype(str) + '_'+ a_clustered.obs[f'rep{i}'].astype(str)
    return adata


from scripts.load_yaml import load
metadata = load('metadata/ggo_config.yml')

adata_dict = dict()
adata_dict2 = dict()

for p in ['PANEL_G', 'PANEL_H']:
    adata_dict[p] = sc.read(metadata[p]['AnnData']['unlabeled_file_name'])
    adata_dict[p] = adata_dict[p][:,~adata_dict[p].var.index.str.contains('EMPTY')]
    adata_dict[p].var['metal'] = adata_dict[p].var.index.str.split('(').str[1].str.replace(')','')
    adata_dict[p].var.index = adata_dict[p].var.index.str.split('(').str[0]
    adata_dict[p].var['panel'] = p

c = 'celltype'
cluster_res = 'cluster_0.5'

# Write Celltype
for panel in adata_dict:

    print(f'Cell type labeling: {panel}...')

    adata = adata_dict[panel]
    adata.obs_names_make_unique()

    adata_dict2[panel] = cluster_rep(adata, n_reps = 2)
    
    # Split celltyps to Epithelial / Stromal / Immune
    adata.obs['celltype_esi'] = 'Etc'
    adata.obs.loc[adata.obs[cluster_res].isin(metadata[f'{panel}_epithelial_stromal']), 'celltype_esi'] = 'Epithelial / Stromal'
    adata.obs.loc[adata.obs[cluster_res].isin(metadata[f'{panel}_immune']), 'celltype_esi'] = 'Immune'
    # adata.write(metadata[f'{panel}_labeled_file_name'])

    # Subset AnnData for more detailed Annotation
    es = adata[adata.obs['celltype_esi'] == 'Epithelial / Stromal']
    immune = adata[adata.obs['celltype_esi'] == 'Immune']

    # es.write(metadata[f'{panel}_epithelial_stromal_name'])
    # immune.write(metadata[f'{panel}_immune_name'])

    adatas = []
    for i, subtype in enumerate(['epithelial_stromal_', 'immune_']):
        labeled_file = metadata[panel]['AnnData'][f'{subtype}labeled_name']
        if not os.path.exists(labeled_file):
            if subtype == 'epithelial_stromal_':
                a = es
            else:
                a = immune
            res = 0.5
            sc.pp.pca(a)
            sc.pp.neighbors(a)
            #sc.tl.umap(a, gamma = 3)
            p = PARC(
                a.X,
                neighbor_graph=a.obsp["connectivities"],
                random_seed=42,
                resolution_parameter=res,
            )
            p.run_PARC()
            a.obs[f"cluster_{res}"] = pd.Categorical(pd.Series(p.labels) + 1)
            a.write(labeled_file)
        else:
            labeled_adata = sc.read(labeled_file)

        celltype_map = metadata[f'{panel}_{subtype}celltype_label_0_5']
        celltype_map2 = {k: f'{k}. {celltype_map[k]}' for k in celltype_map}

        labeled_adata.obs[c] = pd.Categorical(labeled_adata.obs[cluster_res].astype(int).replace(celltype_map))
        labeled_adata.obs['celltype_cid'] = pd.Categorical(labeled_adata.obs[cluster_res].astype(int).replace(celltype_map2))
        labeled_adata = labeled_adata[labeled_adata.obs[c].isin(celltype_map.values()),:]
        labeled_adata.write(labeled_file)
        adatas.append(labeled_adata)


        imc.pl.celltype_heatmap(
            labeled_adata,
            cluster_ids = [cluster_res, 'celltype_cid', 'celltype'],
            var_names = var_names, # metadata[f'{panel}_var_celltype_groups']
            panel = panel + subtype
        )

        umap_adata_filename = metadata[f'{panel}{subtype}umap_name']
        if not os.path.exists(umap_adata_filename):
            umap_adata = sc.pp.subsample(labeled_adata, n_obs = 30000, copy = True)
            sc.pp.pca(umap_adata)
            sc.pp.neighbors(umap_adata)
            sc.tl.umap(umap_adata, gamma = 3)
            umap_adata.write(umap_adata_filename)
        else:
            umap_adata = sc.read(umap_adata_filename)
        sc.pl.umap(umap_adata, color = 'celltype', frameon=False, show=False, save = f'figures/celltype/{panel}{subtype}UMAP_celltype')
        sc.pl.umap(umap_adata, color = 'celltype_cid', frameon=False, show=False, save = f'figures/celltype/{panel}{subtype}UMAP_celltype')

    adata = anndata.concat(adatas)
    adata.obs['celltype_broad'] = adata.obs['celltype'].str.split(' (', regex = False).str[0]
    adata.uns['celltype_colors'] = metadata[f'{panel}_celltype_colors']
    adata.uns['celltype_broad_colors'] = metadata[f'{panel}_celltype_broad_colors']

    adata.obs['Radiology'] = adata.obs['Radiology'].replace('part solid', 'Part Solid').replace('Part solid', 'Part Solid').replace('Can not Identify', 'UNK')
    adata.obs['Radiology'] = pd.Categorical(adata.obs['Radiology'], categories = ['UNK', 'Normal', 'Pure Non Solid', 'Part Solid', 'Solid'])
    adata.obs.loc[adata.obs['sample_y'].str.contains('N'), 'Radiology'] = 'Normal'

    adata.obs['pred'] = adata.obs['pred'].str.upper().replace('LEP/MIA','LEP').replace('LPA', 'LEP').replace('MIA ', 'MIA').replace('MP','Other').replace('SOL','Other')
    adata.obs.loc[adata.obs['sample_y'].str.contains('N'), 'pred'] = 'Normal'
    adata.obs['pred'] = pd.Categorical(adata.obs['pred'], categories = ['Normal', 'AIS', 'MIA', 'LEP', 'AC', 'Other'])
    
    adata.obs['pathology'] = adata.obs['pred'].astype(str)
    adata.obs['pathology'] = adata.obs['pathology'].replace({'AC': 'IAC', 'LEP': 'IAC', 'Other': 'IAC'})
    adata.obs['pathology'] = pd.Categorical(adata.obs['pathology'], categories = ['Normal', 'AIS', 'MIA', 'IAC'])

    adata.obs['radio'] = adata.obs['Radiology'].astype(str)
    adata.obs['radio'] = adata.obs['radio'].replace({'Normal': 'N', 'Pure Non Solid': 'PNS', 'Part Solid': 'PS', 'Solid': 'S'})
    adata.obs['radio'] = pd.Categorical(adata.obs['radio'], categories = ['N', 'PNS', 'PS', 'S', 'UNK'])
    
    df = pd.read_csv('metadata/Panel_G_ROI_area.csv', index_col = 0)
    tmp = adata.obs.merge(df, how = 'left', left_on = 'roi', right_on = 'roi')
    adata.obs['ROI_area'] = tmp['ROI_area'].tolist()

    adata.write(metadata[f'{panel}_phenotyped_file_name'])

    umap_adata_filename = metadata[f'{panel}_phenotyped_umap_name']
    if not os.path.exists(umap_adata_filename):
        umap_adata = sc.pp.subsample(adata, n_obs = 30000, copy = True)
        sc.pp.pca(umap_adata)
        sc.pp.neighbors(umap_adata)
        sc.tl.umap(umap_adata, gamma = 3)

        umap_adata.write(umap_adata_filename)
    else:
        umap_adata = sc.read(umap_adata_filename)
    sc.pl.umap(umap_adata, color = 'celltype', frameon=False, show=False, save = f'figures/celltype/{panel}_UMAP_celltype')
    sc.pl.umap(umap_adata, color = 'celltype_broad', frameon=False, show=False, save = f'figures/celltype/{panel}_UMAP_celltype_broad')


    adata.obs[['GGO ID', 'Radiology', 'pred']].drop_duplicates().groupby(['Radiology', 'pred']).count().reset_index().pivot(index = 'Radiology', columns = 'pred').to_csv(f'metadata/{panel}_radiology_pathology_roi_count.csv')
    adata.obs[['sample_y', 'Radiology', 'pred']].drop_duplicates().groupby(['Radiology', 'pred']).count().reset_index().pivot(index = 'Radiology', columns = 'pred').to_csv(f'metadata/{panel}_radiology_pathology_roi_count.csv')

