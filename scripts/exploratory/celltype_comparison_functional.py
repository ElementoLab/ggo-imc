
import os
from glob import glob
from pathlib import Path
import tifffile
from tqdm import tqdm

import pandas as pd
import numpy as np

import scanpy as sc
import anndata

from skimage.exposure import adjust_gamma
from skimage import filters
import scipy.ndimage as ndi
import scipy

import seaborn as sns
import matplotlib.pyplot as plt

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

c = 'celltype'
cluster_res = 'cluster_0.5'

metadata_filename = 'metadata/ggo_config.yml'

with open(metadata_filename, "r") as stream:
    try:
        metadata = yaml.safe_load(stream)
    except yaml.YAMLError as exc:
        print(exc)

adata_dict = dict()
for p in ['PANEL_G', 'PANEL_H']:
    print(f"Reading {metadata[p]['AnnData']['phenotyped_file_name']}...")
    adata_dict[p] = sc.read(metadata[p]['AnnData']['phenotyped_file_name'])


for panel in ['PANEL_G', 'PANEL_H']:
    for id_ in ['GGO289N', 'GGO289T_2', 'GGO287T_2']:
        fig, ax = plt.subplots(1,1,figsize = (10,3), dpi = 300)
        sc.pl.spatial(
            adata_dict[panel][adata_dict[panel].obs['GGO ID'] == id_],
            color = 'celltype',
            spot_size = 10,
            frameon = False,
            show = False,
            ax = ax
        )
        os.makedirs('figures/spatial/', exist_ok = True)
        plt.tight_layout()
        plt.savefig(f'figures/spatial/{id_}_{panel}_celltype.pdf', bbox_inches = 'tight')
        plt.close()


filt = [{'pathology': 'Normal', 'radio': 'N'}, {'pathology': 'AIS', 'radio': 'PNS'}, {'pathology': 'MIA', 'radio': 'PS'}, {'pathology': 'IAC', 'radio': 'S'}]
for p in adata_dict:
    adata_dict[p] = anndata.concat([
        adata_dict[p][(adata_dict[p].obs['pathology'] == s['pathology']) & (adata_dict[p].obs['radio'] == s['radio'])] for s in filt
    ])
    adata_dict[p] = adata_dict[p][:,~adata_dict[p].var.index.str.contains('EMPTY')]
    adata_dict[p].obs['pathology'] = pd.Categorical(adata_dict[p].obs['pathology'], categories = ['Normal', 'AIS', 'MIA', 'IAC'])
    adata_dict[p].obs['radio'] = pd.Categorical(adata_dict[p].obs['radio'], categories = ['N', 'PNS', 'PS', 'S'])
    adata_dict[p].var.index = adata_dict[p].var.index.str.split('(').str[0]

for p in adata_dict:
    for celltype in metadata['CELLTYPES']:
        for ct in tqdm(adata_dict[p].obs[celltype].unique()):

            a = adata_dict[p][adata_dict[p].obs[celltype] == ct].copy()
            sc.tl.rank_genes_groups(a, groupby = 'radio', method = 'wilcoxon', use_raw = False)
            #sc.pl.rank_genes_groups_dotplot(a, n_genes = 5, values_to_plot = 'logfoldchanges', min_logfoldchange = 2,cmap = 'bwr', show = False, title = celltype, vmax = 5, vmin = -5, dendrogram = False)
            sc.pl.rank_genes_groups_dotplot(a, n_genes = 5, values_to_plot = 'scores', cmap = 'bwr', show = False, title = ct, dendrogram = False, use_raw = False)
            
            path = f'figures/{celltype}_diff/'

            os.makedirs(path, exist_ok = True)
            cts = ct.replace('/','')
            plt.savefig(path + f'dotplot_{p}_{cts}.pdf', bbox_inches = 'tight')
            plt.close()


utag_dict = dict()

for p in ['PANEL_G', 'PANEL_H']:
    print(f"Reading {metadata[p]['AnnData']['utag_labeled_name']}...")
    utag_dict[p] = sc.read(metadata[p]['AnnData']['utag_labeled_name'])

for panel in ['PANEL_G', 'PANEL_H']:
    for id_ in ['GGO289N', 'GGO289T_2', 'GGO287T_2']:
        fig, ax = plt.subplots(1,1,figsize = (10,3), dpi = 300)
        sc.pl.spatial(
            utag_dict[panel][utag_dict[panel].obs['GGO ID'] == id_],
            color = 'uE',
            spot_size = 10,
            frameon = False,
            show = False,
            ax = ax
        )
        os.makedirs('figures/spatial/', exist_ok = True)
        plt.tight_layout()
        plt.savefig(f'figures/spatial/{id_}_{panel}_ue.pdf', bbox_inches = 'tight')
        plt.close()


filt = [{'pathology': 'Normal', 'radio': 'N'}, {'pathology': 'AIS', 'radio': 'PNS'}, {'pathology': 'MIA', 'radio': 'PS'}, {'pathology': 'IAC', 'radio': 'S'}]
for p in utag_dict:
    utag_dict[p] = anndata.concat([
        utag_dict[p][(utag_dict[p].obs['pathology'] == s['pathology']) & (utag_dict[p].obs['radio'] == s['radio'])] for s in filt
    ])
    utag_dict[p] = utag_dict[p][:,~utag_dict[p].var.index.str.contains('EMPTY')]
    utag_dict[p].obs['pathology'] = pd.Categorical(utag_dict[p].obs['pathology'], categories = ['Normal', 'AIS', 'MIA', 'IAC'])
    utag_dict[p].obs['radio'] = pd.Categorical(utag_dict[p].obs['radio'], categories = ['N', 'PNS', 'PS', 'S'])
    utag_dict[p].var.index = utag_dict[p].var.index.str.split('(').str[0]


# for p in utag_dict:
#     adata = utag_dict[p]
#     adata.obs['PID'] = adata.obs['description'].str[:6]

#     for celltype in ['uE', 'uE_broad']:
        
#         density = compute_density(
#             adata,
#             celltype_key = celltype,
#             condition_keys = metadata['CONDITIONS'])
        
#         for cond in metadata['CONDITIONS']:
#             for pval_form in ['star', 'sci_notation']:
#                 if f'{cond}_color' in metadata:
#                     palette = metadata[f'{cond}_color']
#                 else:
#                     palette = None #'tab10'

#                 plot_grouped_key_mwu(
#                     density,
#                     condition_keys = [cond],
#                     palette = palette,
#                     pval_form = pval_form,
#                     save_dir = f'figures/{p}_{celltype}_density_across_condition_filtered/'
#                 )

'09282022_CTMA123F_PanelG_1-09'

for p in utag_dict:
    for celltype in ['uE', 'uE_broad']:
        for ct in tqdm(utag_dict[p].obs[celltype].unique()):

            a = utag_dict[p][utag_dict[p].obs[celltype] == ct].copy()
            sc.tl.rank_genes_groups(a, groupby = 'radio', method = 'wilcoxon', use_raw = False)
            #sc.pl.rank_genes_groups_dotplot(a, n_genes = 5, values_to_plot = 'logfoldchanges', min_logfoldchange = 2,cmap = 'bwr', show = False, title = celltype, vmax = 5, vmin = -5, dendrogram = False)
            sc.pl.rank_genes_groups_dotplot(a, n_genes = 5, values_to_plot = 'scores', cmap = 'bwr', show = False, title = ct, dendrogram = False, use_raw = False)
            
            path = f'figures/{celltype}/'

            os.makedirs(path, exist_ok = True)
            cts = ct.replace('/','')
            plt.savefig(path + f'dotplot_{p}_{cts}.pdf', bbox_inches = 'tight')
            plt.close()

