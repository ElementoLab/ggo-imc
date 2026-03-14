"""Figure 5: patient risk group process scores and conditional probability plots."""

import pandas as pd
import numpy as np
import scanpy as sc
import seaborn as sns
import matplotlib.pyplot as plt
import imc_analysis as imc
from scipy.stats import zscore

from utils import load_config, save_figure, ensure_dir


def cond_prob(data, y, x):
    """Compute the conditional probability P(y | x) as a pivot table.

    Parameters
    ----------
    data:
        DataFrame with columns ``x`` and ``y``.
    y:
        Column whose conditional distribution is computed.
    x:
        Conditioning column (or list of columns).

    Returns
    -------
    pd.DataFrame
        Pivot table with ``x`` as the index and ``y`` categories as columns,
        values are percentages (0–100).
    """
    return (
        data
        .groupby(x)[y]
        .value_counts(normalize=True)
        .mul(100)
        .rename('percent')
        .reset_index()
        .pivot_table(index=x, columns=y, values='percent')
        .reset_index()
    )


if __name__ == '__main__':
    metadata = load_config()

    ensure_dir('figures/figure5')

    adata = sc.read(
        metadata['patient_celltype_broad_clustered'],
        backup_url=metadata['patient_group_url'],
    )

    sc.pp.scale(adata)
    sc.pp.pca(adata)
    sc.pp.neighbors(adata)
    sc.tl.umap(adata)

    group = 'Group'
    features = metadata['CONDITIONS']

    # ------------------------------------------------------------------
    # Process scores from cellular abundance
    # ------------------------------------------------------------------
    processes = {
        'Epithelial score': ['Epi.', 'Epi. Prol.', 'Epi.-like', 'Tumor-like'],
        'Mesenchymal score': ['Mesen.-like'],
        'Fibrosis score': ['Fib.'],
        'Angiogenesis score': ['Endo.'],
        'B&T score': ['B', 'CD4 T', 'CD8 T'],
        'Pan-Immune score': ['CD4 T', 'CD8 T', 'Mac.', 'NK', 'T reg', 'Mast', 'Mono.', 'Neut.', 'PMN-MDSC'],
    }

    adata.obs[group] = adata.obs[group].str.replace('roup', '')
    for s in processes:
        sc.tl.score_genes(adata, gene_list=processes[s], score_name=s, ctrl_size=20, n_bins=10)

    adata.obs['EMT score'] = adata.obs['Mesenchymal score'] - adata.obs['Epithelial score']
    scores = [
        'EMT score', 'Epithelial score', 'Mesenchymal score',
        'Fibrosis score', 'Angiogenesis score', 'B&T score', 'Pan-Immune score',
    ]
    adata.uns['Group_colors'] = ['#99B898', '#FECEAB', '#E84A5F', '#2A363B']

    fig, axes = plt.subplots(1, len(scores), figsize=(12, 2.5), dpi=300)
    for i, ax in enumerate(axes.flatten()):
        adata.obs[scores[i]] = zscore(adata.obs[scores[i]])
        sc.pl.violin(adata, scores[i], groupby=group, stripplot=False, inner='box', show=False, ax=ax)
        ax.set_xlabel('')
    sns.despine()
    plt.tight_layout()
    save_figure('figures/figure5/group_process_scores.pdf')

    # ------------------------------------------------------------------
    # Conditional probability: P(histology | radiology [, group])
    # ------------------------------------------------------------------
    p_hist_radio = cond_prob(adata.obs, x='radio', y='pathology')
    p_hist_radio_group = cond_prob(adata.obs, x=['radio', group], y='pathology')
    groups = p_hist_radio_group.groupby(group)

    fig, axes = plt.subplots(1, len(groups) + 1, sharey=True, figsize=(7, 2.5), dpi=300)

    p_hist_radio.plot(x='radio', kind='bar', stacked=True, ax=axes[0], color=metadata['pathology_color'][1:])
    axes[0].set_xlabel('All Groups')
    axes[0].set_ylabel(r'P($Histology|Radiology$)')
    axes[0].legend().remove()

    for ax, (grp, data) in zip(axes[1:], groups):
        data.set_index('radio').plot(kind='bar', stacked=True, ax=ax, color=metadata['pathology_color'][1:])
        ax.tick_params(axis='both', which='both', length=0)
        ax.set_xlabel(f'{grp}')
        ax.set_ylabel(r'P($Histology|Radiology$)')
        ax.legend().remove()

    sns.despine()
    plt.tight_layout()
    save_figure('figures/figure5/group_conditional_probability_hist_radio.pdf')

    # Information gained: I(histology | radiology, group)
    # = H(histology | radiology) - H(histology | radiology, group)
    print(p_hist_radio)
    print(p_hist_radio_group)
