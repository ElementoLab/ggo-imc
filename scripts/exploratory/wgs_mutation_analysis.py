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

# GGO Variant
df = pd.read_excel('metadata/GGO_variants_Drivers_06022023(1).xlsx')
df['Sample ID'] = df['Tumor_Sample_Barcode'].str[:6]
df.to_csv('metadata/ggo_variant.csv')

# Pivot Table
pdf = df.groupby(['Hugo_Symbol', 'Sample ID']).count()['Chromosome'].reset_index().pivot(index = 'Sample ID', columns = 'Hugo_Symbol')
pdf[pdf > 0] = True
pdf = pdf.fillna(False)
pdf.to_csv('metadata/mutation.csv')
# pdf.to_csv('mutation_piv.csv')

# Patient
adata = sc.read(metadata['PANEL_H']['AnnData']['patient_celltype_name'])
adata.obs['Sample ID'] = adata.obs.index.astype(str).str[:6]

pdf = pd.read_csv('metadata/mutation_piv.csv', index_col = 0)
genes_ = ['EGFR', 'KRAS', 'RBM10', 'TP53']
pdf = pdf[genes_]


obs = adata.obs.merge(pdf, left_on = 'Sample ID', right_index = True, how="left")
idx = adata.obs.index
adata.obs = obs
adata.obs.index = idx

# obs = adata.obs.merge(pdf, left_index = True, right_on = 'Sample ID')
# idx = adata.obs.index.intersection(pdf['Sample ID'])
# adata = adata[idx,:]
# adata.obs = obs
# adata.obs.index = idx


#adata.write('patient_cell_mutation.h5ad')

# adata = sc.read('patient_cell_mutation.h5ad')
adata.var.index = adata.var['group'].astype(str) + '_' + adata.var.index


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



tumor_celltype_list = ['T_B', 'T_CD4 T', 'T_CD8 T', 'T_Endo. (CD31+, aSMA+)', 'T_Epi. (RAGE+)', 'T_Epi.-like (PanCK+)', 'T_Epi.-like (PanCK++)', 'T_Epi.-like (PanCK++, SFTPC+)', 'T_Epi.-like (PanCK+, RAGE+, SFTPC+)', 'T_Fib. (aSMA+)', 'T_Low Expr.', 'T_Mac. (CD68+, CD163-)', 'T_Mac. (CD163+)', 'T_Mast', 'T_Mono.', 'T_Neut.', 'T_PMN-MDSC']
plot_grouped_key_mwu(adata[:,tumor_celltype_list], condition_keys = genes_, save_dir = f'figures/mutation_celltype/tumor/', palette = None)
normal_celltype_list = ['N_B', 'N_CD4 T', 'N_CD8 T', 'N_Endo. (CD31+, aSMA+)', 'N_Epi. (RAGE+)', 'N_Epi.-like (PanCK+)', 'N_Epi.-like (PanCK++)', 'N_Epi.-like (PanCK++, SFTPC+)', 'N_Epi.-like (PanCK+, RAGE+, SFTPC+)', 'N_Fib. (aSMA+)', 'N_Low Expr.', 'N_Mac. (CD68+, CD163-)', 'N_Mac. (CD163+)', 'N_Mast', 'N_Mono.', 'N_Neut.', 'N_PMN-MDSC']
plot_grouped_key_mwu(adata[:,normal_celltype_list], condition_keys = genes_, save_dir = f'figures/mutation_celltype/normal/', palette = None)


adata