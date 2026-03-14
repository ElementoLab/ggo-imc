"""Figure 4: UTAG microenvironment composition, area-based comparison, and marker dotplots."""

import anndata
import scanpy as sc
import seaborn as sns
import matplotlib.pyplot as plt
import imc_analysis as imc
from tqdm import tqdm

from utils import load_config, load_panels, save_figure, ensure_dir


if __name__ == '__main__':
    metadata = load_config()

    adata_dict = load_panels(
        metadata,
        anndata_key='utag_labeled_name',
        backup_key='utag_url',
    )

    ensure_dir('figures/figure4/')

    # ------------------------------------------------------------------
    # Celltype composition per microenvironment (bar charts)
    # ------------------------------------------------------------------
    celltype = 'celltype_broad'
    uE = 'uE_broad'

    for p in ['PANEL_G', 'PANEL_H']:
        count = adata_dict[p].obs.groupby([celltype, uE]).count()
        count = count / 1000

        g = sns.FacetGrid(
            data=count['sample'].reset_index(),
            hue=celltype,
            col=uE,
            col_wrap=2,
            aspect=1.3,
            palette=adata_dict[p].uns['celltype_broad_colors'],
            sharex=False,
            height=3,
        )
        g.map(sns.barplot, 'sample', celltype)

        for ax in g.axes.flat:
            title = ax.get_title().replace(f'{uE} = ', '')
            ax.set_title(title, fontsize=12)
            ax.set_ylabel('', fontsize=12)
            ax.set_xlabel('')

        plt.tight_layout()
        plt.xlabel('x1000 cells')
        save_figure(f'figures/figure4/{p}_{uE}_{celltype}_composition.pdf')

    # ------------------------------------------------------------------
    # Area-based microenvironment comparison (PANEL_G only)
    # ------------------------------------------------------------------
    p = 'PANEL_G'
    density = imc.tl.celltype_density(
        adata_dict[p],
        celltype=uE,
        condition_keys=['pathology', 'radio'],
    )

    area = adata_dict[p].obs.groupby(['roi', uE])[['area']].sum().pivot_table(index='roi', columns=uE)
    area = area / 1e6
    area.columns = area.columns.droplevel(0)

    area_adata = anndata.AnnData(X=area, obs=density.obs, var=density.var)
    imc.tl.grouped_mwu_test(area_adata, condition_keys=['radio', 'pathology'])

    area_adata.uns['radio_colors'] = adata_dict[p].uns['radio_colors']
    area_adata.uns['pathology_colors'] = adata_dict[p].uns['pathology_colors']

    for pval_form in ['star', 'sci_notation']:
        imc.pl.plot_mwu(
            area_adata,
            save_dir=f'figures/figure4/{p}/{uE}_area/',
            pval_form=pval_form,
            kind='violin',
        )

    # ------------------------------------------------------------------
    # Functional marker dotplots per microenvironment and condition
    # ------------------------------------------------------------------
    functional_markers = [
        'HLADR', 'GranzymeB', 'Tbet', 'CD163',
        'ICOS', 'CCR4', 'LAG3', 'CD4',
        'CCR7', 'CD28', 'CD11b',
        'CD45RA', 'CD206', 'Ki67', 'CD57', 'CD27', '41BB', 'CD45RO',
        'TIM3', 'CD25', 'CTLA4',
    ]

    celltype = 'uE_broad'
    for p in ['PANEL_G', 'PANEL_H']:
        for feat in ['pathology', 'radio']:
            for ct in tqdm(adata_dict[p].obs[celltype].unique()):
                panel_markers = [m for m in functional_markers if m in adata_dict[p].var_names]
                a = adata_dict[p][adata_dict[p].obs[celltype] == ct, panel_markers].copy()
                sc.tl.rank_genes_groups(a, groupby=feat, method='wilcoxon', use_raw=False)
                sc.pl.rank_genes_groups_dotplot(
                    a, n_genes=4, values_to_plot='scores',
                    cmap='bwr', show=False, title=ct,
                    dendrogram=False, use_raw=False,
                )

                cts = ct.replace('/', '')
                path = f'figures/figure4/{celltype}/{feat}/'
                ensure_dir(path)
                save_figure(path + f'dotplot_score_{p}_{feat}_{cts}.pdf')
