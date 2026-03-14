"""Figure 3: epithelial phenotype proportions and EMT (PanCK/Vimentin) density plots."""

import pandas as pd
import numpy as np

import scanpy as sc
import seaborn as sns
import matplotlib.pyplot as plt
import imc_analysis as imc
from tqdm import tqdm

from utils import load_config, load_single_panel, save_figure, ensure_dir


if __name__ == '__main__':
    metadata = load_config()

    ensure_dir('figures/figure3/')

    # ------------------------------------------------------------------
    # PANEL_H — epithelial subtype cell proportion
    # ------------------------------------------------------------------
    sc.logging.info(f"Reading {metadata['PANEL_H']['AnnData']['phenotyped_umap_name']}...")
    adata = load_single_panel(metadata, 'PANEL_H')

    celltype_map = {
        'Tumor-like (RAGE+)': 'RAGE+',
        'Tumor-like (SFTPC+)': 'SFTPC+',
        'Tumor-like': 'PanCK+',
        'Epi.-like (RAGE+)': 'RAGE+',
        'Epi.-like (SFTPC+)': 'SFTPC+',
    }

    adata = adata[adata.obs['celltype'].isin(celltype_map.keys())]
    adata.obs['phenotype'] = adata.obs['celltype'].replace(celltype_map)
    adata.obs['phenotype'] = pd.Categorical(
        adata.obs['phenotype'], categories=['PanCK+', 'RAGE+', 'SFTPC+'], ordered=True
    )

    for cond in ['pathology', 'radio']:
        phenotype = adata.obs.groupby(['phenotype', cond]).count()['sample'].unstack()
        phenotype = phenotype / np.array(phenotype.sum(axis=0))[None, :]

        fig, ax = plt.subplots(1, 1, figsize=(3, 2), dpi=300)
        phenotype.T.reset_index().plot(x=cond, kind='bar', stacked=True, title='', ax=ax)

        df_rel = phenotype.T
        df_pos = df_rel.cumsum(1) - df_rel / 2

        for i, x in enumerate(df_rel.index):
            for y in df_rel.columns:
                plt.text(
                    i, df_pos.loc[x, y],
                    str(np.round(df_rel.loc[x, y] * 100, 1)) + '%',
                    va='center', ha='center',
                )

        sns.despine()
        ax.legend().remove()
        ax.set_xlabel('')
        ax.set_xticklabels(ax.get_xticklabels(), rotation=0)
        plt.tight_layout()
        save_figure(f'figures/figure3/epithelial_proportion_{cond}.pdf')

    # ------------------------------------------------------------------
    # PANEL_G — EMT (PanCK vs Vimentin) density plots
    # ------------------------------------------------------------------
    sc.logging.info(f"Reading {metadata['PANEL_G']['AnnData']['phenotyped_umap_name']}...")
    pg = load_single_panel(metadata, 'PANEL_G')

    pg_epi = pg[pg.obs['celltype_broad'].isin(
        ['Epithelial-like', 'Epithelial-like (Ki67+)', 'Mesenchymal-like']
    )]

    fig, axes = plt.subplots(2, 4, dpi=300, figsize=(8, 4))

    for i, ax in tqdm(enumerate(axes.flatten())):
        feature = 'pathology' if int(i / 4) == 0 else 'radio'
        rad = pg_epi.obs[feature].cat.categories[i % 4]
        pg_tmp = pg_epi[pg_epi.obs[feature] == rad].copy()
        sns.kdeplot(pg_tmp.to_df(), x='PanCK', y='Vimentin', fill=True, ax=ax, cmap='YlOrRd')
        ax.set_title(rad)

    plt.tight_layout()
    save_figure('figures/figure3/EMT proportion.pdf')
