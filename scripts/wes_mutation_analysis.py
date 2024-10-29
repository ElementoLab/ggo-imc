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

metadata = load('metadata/ggo_config.yml')

# Patient
adata = sc.read(metadata['PANEL_H']['AnnData']['patient_celltype_name'])
adata.obs['Sample ID'] = adata.obs.index.astype(str).str[:6]

pdf = pd.read_csv(metadata['wes_variant_pivot_saved'], index_col = 0)
genes_ = ['EGFR', 'KRAS', 'RBM10', 'TP53']
pdf = pdf[genes_]

# concatenate exome sequencing mutation
obs = adata.obs.merge(pdf, left_on = 'Sample ID', right_index = True, how="left")
idx = adata.obs.index
adata.obs = obs
adata.obs.index = idx



# Mutation proportion
gene = ['ARID1A', 'ATM', 'BRAF', 'CDKN2A', 'CTNNB1', 'EGFR', 'ERBB2', 'KEAP1', 'KRAS', 'MET', 'MGA', 'NF1', 'PIK3CA', 'RB1', 'RBM10', 'RIT1', 'SETD2', 'SMARCA4', 'STK11', 'TP53', 'U2AF1']
print(adata.obs[gene].astype(float).sum(axis = 0).sort_values(ascending = False))
for g in gene:
    adata.obs[g] = pd.Categorical(adata.obs[g])

mut = adata.obs[gene].astype(float).sum(axis = 0).sort_values().tail(4).to_frame().rename(columns={0: 'MUT'})#.reset_index()
mut['WT'] = 123-17-mut['MUT']
mut['NA'] = 17.0

fig, ax = plt.subplots(1,1,dpi = 300, figsize = (3,2))
mut.plot.barh(legend = False, ax = ax, color = ['#DCCDE8','#EFC88B',  '#BBBBBB'], stacked = True)
for container in ax.containers:
    ax.bar_label(container)
sns.despine()
ax.set_xlabel('Count')
ax.set_ylabel('Gene')
ax.set_title('')
#ax.set_xlim(0,50)
plt.tight_layout()
plt.savefig('figures/mutation_count_color.pdf')
plt.close()
#plt.show()

