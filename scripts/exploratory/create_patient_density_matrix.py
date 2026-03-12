
import os
from tqdm import tqdm

import pandas as pd
import numpy as np

import scanpy as sc
import anndata

import seaborn as sns
import matplotlib.pyplot as plt

import anndata

import imc_analysis as imc
from scripts.load_yaml import load


# load data
metadata = imc.utils.parse_yaml('metadata/ggo_config.yml')

pdf = pd.read_csv(metadata['wes_variant_pivot_saved'], index_col = 0)
genes_ = ['EGFR', 'KRAS', 'RBM10', 'TP53']
pdf = pdf[genes_]

# Read in data
adata_dict = dict()
for p in ['PANEL_G', 'PANEL_H']:
    print(f"Reading {metadata[p]['AnnData']['phenotyped_file_name']}...")
    adata_dict[p] = sc.read(metadata[p]['AnnData']['phenotyped_file_name'])
    adata_dict[p].obs['PID'] = adata_dict[p].obs['description'].str[:6]
    adata_dict[p].obs['radio'] = pd.Categorical(adata_dict[p].obs['radio'], categories = ['N', 'PNS', 'PS', 'S'], ordered=True)
    adata_dict[p] = adata_dict[p][adata_dict[p].obs['radio'] != 'N']
    adata_dict[p].obs['pathology'] = pd.Categorical(adata_dict[p].obs['pathology'], categories = ['Normal', 'AIS', 'MIA', 'IAC'], ordered = True)
    del adata_dict[p].obs['ROI_area']

    # attach wes
    obs = adata_dict[p].obs.merge(pdf, left_on = 'PID', right_index = True, how="left")
    idx = adata_dict[p].obs.index
    adata_dict[p].obs = obs
    adata_dict[p].obs.index = idx

    
    for g in genes_:
        adata_dict[p].obs[g] = pd.Categorical(adata_dict[p].obs[g].replace({False: 'WT', True:'MUT'}), categories = ['WT', 'MUT'])

# Get cell type density matrix
celltype = 'celltype_broad'
p_density = dict()

for p in ['PANEL_G', 'PANEL_H']:
    adata = adata_dict[p].copy()
    p_density[p] = imc.tl.patient_density(
        adata,
        celltype_key = celltype,
        patient_key = 'PID',
        condition_keys = metadata['CONDITIONS'] + genes_ + ['Smoking Status', 'Gender', 'Race'])

    if 'GGO ID' in p_density[p].obs:
        idx_keep = p_density[p].obs[['GGO ID']].drop_duplicates().index
        p_density[p] = p_density[p][idx_keep]

        p_density[p].obs.set_index('GGO ID', inplace = True)

# Get common rois and celltypes
common_roi = pd.Index.intersection(*[p_density[p].obs.index for p in p_density])
celltype_drop = ['Low Expr.', 'Other', 'B / Memory T', 'Naive B']

# Save panels
for p in p_density:
    p_density[p] = p_density[p][common_roi,~p_density[p].var.index.isin(celltype_drop)]

common_celltype = pd.Index.intersection(*[p_density[p].var.index for p in p_density])
for ct in common_celltype:
    p_density['PANEL_G'][:,ct].X = p_density['PANEL_G'][:,ct].X + p_density['PANEL_H'][:,ct].X / 2
    p_density['PANEL_H'] = p_density['PANEL_H'][:,p_density['PANEL_H'].var.index != ct]

# concatenate densities
density = anndata.concat([p_density[p].copy() for p in p_density], axis = 1)
density.obs = p_density['PANEL_G'].obs
density.obs['ROI_area'] = density.obs['ROI_area'] + p_density['PANEL_H'].obs['ROI_area']

density.write(metadata['patient_celltype_broad_density'])