import pandas as pd
import numpy as np

import scanpy as sc

import seaborn as sns
import matplotlib.pyplot as plt

import imc_analysis as imc

metadata = imc.utils.parse_yaml('metadata/ggo_config.yml')

print(f"Reading {metadata['PANEL_H']['AnnData']['phenotyped_umap_name']}...")
adata = sc.read(
    metadata['PANEL_H']['AnnData']['phenotyped_umap_name'],
    backup_url = metadata['PANEL_H']['AnnData']['backup_url'])

celltype_map = {
    'Tumor-like (RAGE+)':'RAGE+',
    'Tumor-like (SFTPC+)':'SFTPC+',
    'Tumor-like':'PanCK+',
    'Epi.-like (RAGE+)':'RAGE+',
    'Epi.-like (SFTPC+)':'SFTPC+',
}

# Cell Proportion Figure
adata = adata[adata.obs['celltype'].isin(celltype_map.keys())]
adata.obs['phenotype'] = adata.obs['celltype'].replace(celltype_map)
adata.obs['phenotype'] = pd.Categorical(adata.obs['phenotype'], categories = ['PanCK+', 'RAGE+', 'SFTPC+'], ordered = True)
for cond in ['pathology', 'radio']:
    phenotype = adata.obs.groupby(['phenotype', cond]).count()['sample'].unstack()
    phenotype = phenotype / np.array(phenotype.sum(axis = 0))[None,:]

    fig, ax = plt.subplots(1,1,figsize = (3,2), dpi = 300)
    phenotype.T.reset_index().plot(
        x=cond,
        kind='bar',
        stacked=True,
        title='',
        ax = ax)

    df_rel = phenotype.T
    df_pos = df_rel.cumsum(1)-df_rel/2

    for i, x in enumerate(df_rel.index):
        for y in df_rel.columns:
            # print(x, df_pos.loc[x,y], str(np.round(df_rel.loc[x,y] *100, 1)))
            plt.text(i, df_pos.loc[x,y], str(np.round(df_rel.loc[x,y] *100, 1)) + '%', va = 'center', ha = 'center')
    sns.despine()
    ax.legend().remove()
    ax.set_xlabel('')
    ax.set_xticklabels(ax.get_xticklabels(), rotation = 0)
    plt.tight_layout()
    plt.savefig(f'figures/epithelial_proportion_{cond}.pdf')
    plt.close()


# EMT Proportion Figure
print(f"Reading {metadata['PANEL_G']['AnnData']['phenotyped_umap_name']}...")
pg = sc.read(
    metadata['PANEL_G']['AnnData']['phenotyped_umap_name'],
    backup_url = metadata['PANEL_G']['AnnData']['backup_url'])

from tqdm import tqdm
pg_epi = pg[pg.obs['celltype_broad'].isin(['Epithelial-like', 'Epithelial-like (Ki67+)','Mesenchymal-like'])]
#fig, axes = plt.subplots(2,5, dpi = 300, figsize = (10,4))
fig, axes = plt.subplots(2,4, dpi = 300, figsize = (8,4))

for i, ax in tqdm(enumerate(axes.flatten())):
    if int(i/4) == 0:
        feature = 'radio' #'Radiology'
    else:
        feature = 'pathology' #'pred'
    rad = pg_epi.obs[feature].cat.categories[i % 4]
    pg_tmp = pg_epi[pg_epi.obs[feature] == rad].copy()
    sns.kdeplot(pg_tmp.to_df(),
    x = 'PanCK', y = 'Vimentin', fill = True, ax = ax, cmap = 'YlOrRd')
    ax.set_title(rad)

plt.tight_layout()
plt.savefig(f'figures/EMT proportion.pdf', bbox_inches = 'tight')
plt.close()



