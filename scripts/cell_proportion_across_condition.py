
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

import matplotlib
sc.settings.set_figure_params(dpi=200, dpi_save=300, fontsize=12)
matplotlib.rcParams["pdf.fonttype"] = 42
matplotlib.rcParams["ps.fonttype"] = 42
matplotlib.rcParams["axes.grid"] = False
matplotlib.use('Agg')


import yaml

def plot_condition_mwu(
    data: anndata.AnnData,
    figsize: tuple = (3.5,3),
    condition_key: str = 'Radiology',
    sample_key: str = 'roi',
    fudge_key: str = 'sample',
    save_dir: str = None, #'figures/ROI cell count ts.pdf',
    pval_form: str = 'star',
    plot_title: str = '',
    line_width: float = 1,
    verbose: bool = True,
    palette = None
) -> pd.DataFrame:
    import pingouin as pg
    from itertools import combinations
    if pval_form != 'star' and pval_form != 'sci_notation':
        print(f"'pval_form' needs to be either 'star' or 'sci_notation'. Defaulting into 'star'.")
        pval_form = 'star'
        
    def pval_to_star(pval):
        if pval < 0.0001:
            return ' **** '
        elif pval < 0.001:
            return ' *** '
        elif pval < 0.01:
            return ' ** '
        elif pval < 0.05:
            return ' * '
        else:
            return ' ns '
    
    def pval_to_sci_not(pval):
        return "{:.2E}".format(pval)
        
    fig, ax = plt.subplots(1,1, figsize = figsize)
    df = data.obs.groupby([condition_key,sample_key]).count()[fudge_key].reset_index()
    df = df[df[fudge_key]!=0]

    if palette:
        sns.boxplot(
            data = df,
            x = condition_key,
            y = fudge_key,
            fliersize = 0,
            palette = palette,
            ax = ax
        )
    else:
        sns.boxplot(
            data = df,
            x = condition_key,
            y = fudge_key,
            fliersize = 0,
            ax = ax
        )
    sns.swarmplot(
        data = df,
        x = condition_key,
        y = fudge_key,
        s = 3,
        color = 'black',
        ax = ax
    )

    res = list(combinations(df[condition_key].cat.categories, 2))

    pvals = pd.concat([pg.mwu(df[df[condition_key] == p1][fudge_key].tolist(), df[df[condition_key] == p2][fudge_key].tolist()) for p1, p2 in res])
    BH_pvals = pg.multicomp(pvals['p-val'].tolist(), method = 'BH')
    BH_pvals = pd.DataFrame({'Significant': BH_pvals[0], 'adj. p-val': BH_pvals[1]}, index = res)
    
    if verbose:
        print(BH_pvals)
    
    y1 = df[fudge_key].max() * 1.05
    r = y1 * 0.02
    l = df[condition_key].cat.categories.tolist()
    
    if len(BH_pvals) == 1:
        p = BH_pvals['adj. p-val'].tolist()[0]
        if pval_form == 'star':
            pval = pval_to_star(p)
        else:
            pval = pval_to_sci_not(p)
        ax.text(s = pval, x = 0.5, y = y1, fontsize = 8, va = 'bottom', ha = 'center')
        
    sig_n = 0
    for i, sig in enumerate(BH_pvals.index):
        
        if len(BH_pvals) == 1:
            p = BH_pvals['adj. p-val'].tolist()[0]
        else:
            p = BH_pvals.iloc[i]['adj. p-val']

        if p < 0.05:
            ax.plot([l.index(sig[0]), l.index(sig[1])], [y1 + r*i, y1 + r*i], lw=line_width, c='black')
            if pval_form == 'star':
                pval = pval_to_star(p)
            else:
                pval = pval_to_sci_not(p)
            ax.text(s = pval, x = 0.5, y = y1, fontsize = 8, va = 'bottom', ha = 'center')

            ax.plot([l.index(sig[0]), l.index(sig[1])], [y1 + r*sig_n, y1 + r*sig_n], lw=line_width, c='black')
            if pval_form == 'star':
                pval = pval_to_star(p)
            else:
                pval = pval_to_sci_not(p)
            ax.text(s = pval, x = l.index(sig[1]), y = y1 + r*sig_n, fontsize = 8, va = 'top', ha = 'left')
            sig_n += 1
    
    ax.grid(visible=False)
    sns.despine()
    ax.set_xlabel('')
    ax.set_ylabel('')
    ax.set_title(plot_title)
    plt.tight_layout()
    if save_dir:
        plt.savefig(save_dir, bbox_inches = 'tight')
        
    
    if verbose:
        plt.show()
    if ax:
        return ax
    plt.close()
    return BH_pvals

# def compute_density(
#     data: anndata.AnnData,
#     condition_key: str = 'Radiology',
#     celltype_key: str = 'celltype',
#     obs_key: str = 'roi',
# ):
    
#     assert(obs_key in data.obs and condition_key in data.obs and celltype_key in data.obs)
#     cell_counts = data.obs.groupby([obs_key,celltype_key]).count()['sample'].reset_index().pivot(index = 'roi', columns = celltype_key, values = 'sample')
#     data_meta = data.obs[[obs_key, condition_key, 'ROI_area']].drop_duplicates()[[obs_key, 'ROI_area', condition_key]]
#     density = cell_counts.merge(data_meta[[obs_key, 'ROI_area']], left_index = True, right_on = obs_key)

#     import numpy as np

#     if obs_key in density:
#         density = density.set_index(obs_key)
        
#     density = density / density['ROI_area'][:,np.newaxis]
#     density

#     if 'ROI_area' in density:
#         del density['ROI_area']
#     density = density.merge(data_meta[[obs_key, condition_key]], left_index = True, right_on = obs_key)
#     return density


def compute_density(
    data: anndata.AnnData,
    condition_keys: list = [],
    celltype_key: str = 'celltype',
    obs_key: str = 'roi',
):
    keys = [obs_key, 'ROI_area'] + condition_keys
    for k in keys:
        assert(k in data.obs)

    cell_counts = data.obs.groupby([obs_key,celltype_key]).count()['sample'].reset_index().pivot(index = obs_key, columns = celltype_key, values = 'sample')
    data_meta = data.obs[keys].drop_duplicates()
    
    # assert cell count matches metadata length
    assert(len(cell_counts) == len(data_meta))
    
    # set index as obs_key
    data_meta = data_meta.set_index(obs_key)
    # cell_counts.index = data_meta.index

    density = anndata.AnnData(X = cell_counts.astype(np.float32), obs = data_meta.loc[cell_counts.index])
    density.var.index = cell_counts.columns

    if obs_key in density.obs:
        density.obs = density.obs.set_index(obs_key)

    density.X = density.X / np.array(density.obs['ROI_area'])[:, None]
    
    return density

def roi_density_to_patient_count(
    adata: anndata.AnnData,
    patient_key: str = 'ID',
    roi_area: str = 'ROI_area'
):
    assert(patient_key in adata.obs)
    assert(roi_area in adata.obs)

    density = adata.to_df()
    count = density * np.array(adata.obs[roi_area])[:,None]
    count[patient_key] = adata.obs[patient_key]
    count = count.groupby(patient_key).sum()
    #count = count.set_index(patient_key)

    features = (adata.obs.groupby(patient_key).nunique() <= 1).all()
    patient_metadata = adata.obs[features[features].index.tolist() + [patient_key]].drop_duplicates()
    patient_metadata = patient_metadata.set_index(patient_key)

    res = anndata.AnnData(X = count, obs = patient_metadata.loc[count.index])
    res.obs[roi_area] = adata.obs.groupby(patient_key)[roi_area].sum().tolist()

    return res

def grouped_obs_mean(adata, group_key, layer=None, gene_symbols=None):
    if layer is not None:
        getX = lambda x: x.layers[layer]
    else:
        getX = lambda x: x.X
    if gene_symbols is not None:
        new_idx = adata.var[idx]
    else:
        new_idx = adata.var_names

    grouped = adata.obs.groupby(group_key)
    out = pd.DataFrame(
        np.zeros((adata.shape[1], len(grouped)), dtype=np.float64),
        columns=list(grouped.groups.keys()),
        index=adata.var_names
    )

    for group, idx in grouped.indices.items():
        X = getX(adata[idx])
        out[group] = X.mean(axis=0, dtype=np.float64).tolist()
    return out.T

def compute_celltype_mean(
    data: anndata.AnnData,
    condition_keys: list = [],
    celltype_key: str = 'celltype',
    obs_key: str = 'roi',
):
    keys = [obs_key, 'ROI_area'] + condition_keys
    for k in keys:
        assert(k in data.obs)

    cell_counts = data.obs.groupby([obs_key,celltype_key]).count()['sample'].reset_index().pivot(index = obs_key, columns = celltype_key, values = 'sample')
    data_meta = data.obs[keys].drop_duplicates()
    
    # assert cell count matches metadata length
    assert(len(cell_counts) == len(data_meta))
    
    # set index as obs_key
    data_meta = data_meta.set_index(obs_key)
    # cell_counts.index = data_meta.index

    density = anndata.AnnData(X = cell_counts.astype(np.float32), obs = data_meta.loc[cell_counts.index])
    density.var.index = cell_counts.columns

    if obs_key in density.obs:
        density.obs = density.obs.set_index(obs_key)

    density.X = density.X / np.array(density.obs['ROI_area'])[:, None]
    
    return density

def compute_mean(
    data: anndata.AnnData,
    condition_keys: list = [],
    obs_key: str = 'roi',
):
    keys = [obs_key] + condition_keys
    for k in keys:
        assert(k in data.obs)

    mean = grouped_obs_mean(data, group_key = obs_key)
    data_meta = data.obs[keys].drop_duplicates()
    
    # assert cell count matches metadata length
    assert(len(mean) == len(data_meta))
    
    # set index as obs_key
    data_meta = data_meta.set_index(obs_key)
    # mean.index = data_meta.index

    mean_adata = anndata.AnnData(X = mean.astype(np.float32), obs = data_meta.loc[mean.index])
    mean_adata.var.index = mean.columns

    if obs_key in mean_adata.obs:
        mean_adata.obs = mean_adata.obs.set_index(obs_key)
    
    return mean_adata

def plot_grouped_key_mwu(
    adata: anndata.AnnData,
    condition_keys: list = None,
    line_width:float = 0.5,
    save_dir: Path = 'figures/celltype_count.pdf',
    palette: list = [],
    pval_form: str = 'star',
    verbose: bool = False
):

    import pingouin as pg
    
    if condition_keys == None:
        condition_keys = adata.obs.columns

    for key in condition_keys:
        assert(key in adata.obs.columns)
        #assert(isinstance(adata.obs[key].dtype, pd.CategoricalDtype))
    
    density = adata.to_df()

    for cond in tqdm(condition_keys):
        
        fig, axes = plt.subplots(3,6, dpi = 300, figsize = (12,8))
        
        for i, ax in enumerate(axes.flatten()):
        
            if i >= density.shape[1]:
                ax.axis('off')
            else:
                
                ct = density.columns[i]

                sns.boxplot(y = density[ct], x = adata.obs[cond], hue = adata.obs[cond], ax = ax, fliersize = 0, palette = palette, dodge = False)
                sns.swarmplot(y = density[ct], x = adata.obs[cond], color = 'black', ax = ax, s = 3)
                ax.set_title(ct)
                ax.set_ylabel('')
                ax.set_xlabel('')

                from itertools import combinations

                def pval_to_star(pval):
                    if pval < 0.0001:
                        return ' **** '
                    elif pval < 0.001:
                        return ' *** '
                    elif pval < 0.01:
                        return ' ** '
                    elif pval < 0.05:
                        return ' * '
                    else:
                        return ' ns '

                def pval_to_sci_not(pval):
                    return "{:.2E}".format(pval)

                res = list(combinations(adata.obs[cond].cat.categories, 2))

                pvals = pd.concat([pg.mwu(density[adata.obs[cond] == p1][ct].tolist(), density[adata.obs[cond] == p2][ct].tolist()) for p1, p2 in res])
                BH_pvals = pg.multicomp(pvals['p-val'].tolist(), method = 'BH')
                BH_pvals = pd.DataFrame({'Significant': BH_pvals[0], 'adj. p-val': BH_pvals[1]}, index = res)

                if verbose:
                    print(BH_pvals)

                y1 = density[ct].max() * 1.05
                r = y1 * 0.03
                l = adata.obs[cond].cat.categories.tolist()

                if len(BH_pvals) == 1:
                    if pval_form == 'star':
                        pval = pval_to_star(BH_pvals['adj. p-val'].tolist()[0])
                    else:
                        pval = pval_to_sci_not(BH_pvals['adj. p-val'].tolist()[0])
                    ax.text(s = pval, x = 0.5, y = y1, fontsize = 8, va = 'bottom', ha = 'center')

                sig_n = 0
                for i, sig in enumerate(BH_pvals.index):

                    if len(BH_pvals) == 1:

                        ax.plot([l.index(sig[0]), l.index(sig[1])], [y1 + r*i, y1 + r*i], lw=line_width, c='black')
                        if pval_form == 'star':
                            pval = pval_to_star(BH_pvals['adj. p-val'].tolist()[0])
                        else:
                            pval = pval_to_sci_not(BH_pvals['adj. p-val'].tolist()[0])
                        ax.text(s = pval, x = 0.5, y = y1, fontsize = 8, va = 'bottom', ha = 'center')

                    else:
                        p = BH_pvals.iloc[i]['adj. p-val']
                        if p < 0.05:
                            ax.plot([l.index(sig[0]), l.index(sig[1])], [y1 + r*sig_n, y1 + r*sig_n], lw=line_width, c='black')
                            if pval_form == 'star':
                                pval = pval_to_star(BH_pvals['adj. p-val'].tolist()[i])
                            else:
                                pval = pval_to_sci_not(BH_pvals['adj. p-val'].tolist()[i])
                            ax.text(s = pval, x = l.index(sig[1]), y = y1 + r*sig_n, fontsize = 8, va = 'top', ha = 'left')
                            sig_n += 1
                ax.legend().remove()
        plt.tight_layout()

        dir_path = f'{save_dir}/{cond}_{pval_form}.pdf'
        # check if directory exists
        if not os.path.exists(save_dir):
            # create directory if it doesn't exist
            os.makedirs(save_dir)
            print(f"Directory '{save_dir}' created.")
        plt.savefig(dir_path, bbox_inches = 'tight')
        if verbose:
            plt.show()
        plt.close()

# def plot_grouped_key_mwu(
#     density: pd.DataFrame,
#     condition_key: str,
#     line_width:float = 0.5,
#     save_key: Path = 'figures/celltype_count.pdf',
#     palette: list = [],
#     pval_form: str = 'star',
#     verbose: bool = False
# ):
#     import pingouin as pg
#     assert(condition_key in density)
#     density = density.reset_index(drop = True)

#     fig, axes = plt.subplots(4,4, dpi = 300, figsize = (12,12))
#     for i, ax in enumerate(axes.flatten()):
#         if i >= len(density.columns) - 1:
#             ax.axis('off')
#         else:
            
#             ct = density.columns[i]

#             sns.boxplot(data = density, y = ct, x = condition_key, ax = ax, fliersize = 0, palette = palette)
#             sns.swarmplot(data = density, y = ct, x = condition_key, color = 'black', ax = ax, s = 3)
#             ax.set_title(ct)
#             ax.set_ylabel('')
#             ax.set_xlabel('')

#             from itertools import combinations

#             def pval_to_star(pval):
#                 if pval < 0.0001:
#                     return ' **** '
#                 elif pval < 0.001:
#                     return ' *** '
#                 elif pval < 0.01:
#                     return ' ** '
#                 elif pval < 0.05:
#                     return ' * '
#                 else:
#                     return ' ns '

#             def pval_to_sci_not(pval):
#                 return "{:.2E}".format(pval)

#             res = list(combinations(density[condition_key].cat.categories, 2))

#             pvals = pd.concat([pg.mwu(density[density[condition_key] == p1][ct].tolist(), density[density[condition_key] == p2][ct].tolist()) for p1, p2 in res])
#             BH_pvals = pg.multicomp(pvals['p-val'].tolist(), method = 'BH')
#             BH_pvals = pd.DataFrame({'Significant': BH_pvals[0], 'adj. p-val': BH_pvals[1]}, index = res)

#             if verbose:
#                 print(BH_pvals)

#             y1 = density[ct].max() * 1.05
#             r = y1 * 0.03
#             l = density[condition_key].cat.categories.tolist()

#             if len(BH_pvals) == 1:
#                 if pval_form == 'star':
#                     pval = pval_to_star(BH_pvals['adj. p-val'].tolist()[0])
#                 else:
#                     pval = pval_to_sci_not(BH_pvals['adj. p-val'].tolist()[0])
#                 ax.text(s = pval, x = 0.5, y = y1, fontsize = 8, va = 'bottom', ha = 'center')

#             sig_n = 0
#             for i, sig in enumerate(BH_pvals.index):

#                 if len(BH_pvals) == 1:

#                     ax.plot([l.index(sig[0]), l.index(sig[1])], [y1 + r*i, y1 + r*i], lw=line_width, c='black')
#                     if pval_form == 'star':
#                         pval = pval_to_star(BH_pvals['adj. p-val'].tolist()[0])
#                     else:
#                         pval = pval_to_sci_not(BH_pvals['adj. p-val'].tolist()[0])
#                     ax.text(s = pval, x = 0.5, y = y1, fontsize = 8, va = 'bottom', ha = 'center')

#                 else:
#                     p = BH_pvals.iloc[i]['adj. p-val']
#                     if p < 0.05:
#                         ax.plot([l.index(sig[0]), l.index(sig[1])], [y1 + r*sig_n, y1 + r*sig_n], lw=line_width, c='black')
#                         if pval_form == 'star':
#                             pval = pval_to_star(BH_pvals['adj. p-val'].tolist()[i])
#                         else:
#                             pval = pval_to_sci_not(BH_pvals['adj. p-val'].tolist()[i])
#                         ax.text(s = pval, x = l.index(sig[1]), y = y1 + r*sig_n, fontsize = 8, va = 'top', ha = 'left')
#                         sig_n += 1
#     plt.tight_layout()

#     plt.savefig(save_key, bbox_inches = 'tight')
#     if verbose:
#         plt.show()
#     plt.close()

# for cond in metadata['CONDITIONS']:
#     if f'{cond}_color' in metadata:
#         palette = metadata[f'{cond}_color']
#     else:
#         palette = None #'tab10'
#     for pval_form in ['star', 'sci_notation']:
#         pvals = plot_condition_mwu(
#             data = pg_phenotyped,
#             condition_key = cond,
#             save_dir = f'figures/ROI cell count {cond} panel g {pval_form}.pdf',
#             figsize = (4.5,3),
#             pval_form = pval_form,
#             palette = palette,
#             verbose = False
#         )

#         pvals = plot_condition_mwu(
#             data = ph_phenotyped,
#             condition_key = cond,
#             save_dir = f'figures/ROI cell count {cond} panel h {pval_form}.pdf',
#             figsize = (4.5,3),
#             pval_form = pval_form,
#             palette = palette,
#             verbose = False
#         )



# for cond in metadata['CONDITIONS']:

#     for celltype in metadata['CELLTYPES']:
#         pg_density = compute_density(pg_phenotyped, celltype_key = celltype, condition_key = cond)
#         ph_density = compute_density(ph_phenotyped, celltype_key = celltype, condition_key = cond)

#         celltypes_g = pg_phenotyped.obs[celltype].cat.categories
#         celltypes_h = ph_phenotyped.obs[celltype].cat.categories

#         epithelial_celltypes_g = celltypes_g[celltypes_g.str.contains('Epi')].tolist()
#         stromal_celltypes_g = celltypes_g[celltypes_g.str.contains('Endo') | celltypes_g.str.contains('Fib') | celltypes_g.str.contains('Mesen')].tolist()
#         immune_celltypes_g = celltypes_g[~(celltypes_g.isin(epithelial_celltypes_g) | celltypes_g.isin(stromal_celltypes_g) | celltypes_g.str.contains('Other'))].tolist()

#         epithelial_celltypes_h = celltypes_h[celltypes_h.str.contains('Epi')].tolist()
#         stromal_celltypes_h = celltypes_h[celltypes_h.str.contains('Endo') | celltypes_h.str.contains('Fib') | celltypes_h.str.contains('Mesen')].tolist()
#         immune_celltypes_h = celltypes_h[~(celltypes_h.isin(epithelial_celltypes_h) | celltypes_h.isin(stromal_celltypes_h) | celltypes_h.str.contains('Other'))].tolist()
        
#         epithelial_g = epithelial_celltypes_g + [cond]
#         stromal_g = stromal_celltypes_g + [cond]
#         immune_g = immune_celltypes_g + [cond]
#         epithelial_h = epithelial_celltypes_h + [cond]
#         stromal_h = stromal_celltypes_h + [cond]
#         immune_h = immune_celltypes_h + [cond]

#         pg_density_epithelial = pg_density[epithelial_g]
#         pg_density_stromal = pg_density[stromal_g]
#         pg_density_immune = pg_density[immune_g]

#         ph_density_epithelial = ph_density[epithelial_h]
#         ph_density_stromal = ph_density[stromal_h]
#         ph_density_immune = ph_density[immune_h]

#         for pval_form in ['star', 'sci_notation']:
#             if f'{cond}_color' in metadata:
#                 palette = metadata[f'{cond}_color']
#             else:
#                 palette = None #'tab10'

#             plot_grouped_key_mwu(
#                 pg_density_epithelial,
#                 condition_key = cond,
#                 palette = palette,
#                 pval_form = pval_form,
#                 save_key = f'figures/cell_density_across_condition/{celltype}_panel_g_{cond}_epithelial_{pval_form}.pdf'
#             )

#             plot_grouped_key_mwu(
#                 ph_density_epithelial,
#                 condition_key = cond,
#                 palette = palette,
#                 pval_form = pval_form,
#                 save_key = f'figures/cell_density_across_condition/{celltype}_panel_h_{cond}_epithelial_{pval_form}.pdf'
#             )
            
#             plot_grouped_key_mwu(
#                 pg_density_stromal,
#                 condition_key = cond,
#                 palette = palette,
#                 pval_form = pval_form,
#                 save_key = f'figures/cell_density_across_condition/{celltype}_panel_g_{cond}_stromal_{pval_form}.pdf'
#             )

#             plot_grouped_key_mwu(
#                 ph_density_stromal,
#                 condition_key = cond,
#                 palette = palette,
#                 pval_form = pval_form,
#                 save_key = f'figures/cell_density_across_condition/{celltype}_panel_h_{cond}_stromal_{pval_form}.pdf'
#             )
            
#             plot_grouped_key_mwu(
#                 pg_density_immune,
#                 condition_key = cond,
#                 palette = palette,
#                 pval_form = pval_form,
#                 save_key = f'figures/cell_density_across_condition/{celltype}_panel_g_{cond}_immune_{pval_form}.pdf'
#             )

#             plot_grouped_key_mwu(
#                 ph_density_immune,
#                 condition_key = cond,
#                 palette = palette,
#                 pval_form = pval_form,
#                 save_key = f'figures/cell_density_across_condition/{celltype}_panel_h_{cond}_immune_{pval_form}.pdf'
#             )
#             

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


conditions = metadata['CONDITIONS'] + ['PID','description']

for p in ['PANEL_G', 'PANEL_H']:
    
    adata = adata_dict[p]
    adata.obs['PID'] = adata.obs['description'].str[:6]

    for celltype in metadata['CELLTYPES']:
        
        density = compute_density(
            adata,
            celltype_key = celltype,
            condition_keys = conditions)
        
        density.write(metadata[p]['AnnData'][f'roi_{celltype}_name'])

        d = {
            'normal_': density[density.obs['radio'] == 'N'],
            'tumor_': density[density.obs['radio'] != 'N'],
        }
        for k in d:
            patient_cellcount = roi_density_to_patient_count(d[k], patient_key = 'PID')
            if 'description' in patient_cellcount.obs:
                del patient_cellcount.obs['description']
            patient_cellcount.write(metadata[p]['AnnData'][f'patient_{k}{celltype}_name'])
            print(p, k, celltype, patient_cellcount.obs.shape)

        for cond in density.obs:
            density.obs[cond] = pd.Categorical(density.obs[cond])

        for pval_form in ['star', 'sci_notation']:
            if f'{cond}_color' in metadata:
                palette = metadata[f'{cond}_color']
            else:
                palette = None #'tab10'

            plot_grouped_key_mwu(
                density,
                condition_keys = metadata['CONDITIONS'],
                palette = palette,
                pval_form = pval_form,
                save_dir = f'figures/{p}_{celltype}_density_across_condition/'
            )

'''
repeat the same, but for samples that are N, PNS-AIS, PS-MIA, S-IAC
'''
for p in ['PANEL_G', 'PANEL_H']:
    
    adata = adata_dict[p]
    adata.obs['PID'] = adata.obs['description'].str[:6]

    for celltype in metadata['CELLTYPES']:
        
        density = compute_density(
            adata,
            celltype_key = celltype,
            condition_keys = conditions)
        
        filt = [{'pathology': 'Normal', 'radio': 'N'}, {'pathology': 'AIS', 'radio': 'PNS'}, {'pathology': 'MIA', 'radio': 'PS'}, {'pathology': 'IAC', 'radio': 'S'}]
        
        density = anndata.concat([
            density[(density.obs['pathology'] == s['pathology']) & (density.obs['radio'] == s['radio'])] for s in filt
        ])
        density.obs['pathology'] = pd.Categorical(density.obs['pathology'], categories = ['Normal', 'AIS', 'MIA', 'IAC'])
        density.obs['radio'] = pd.Categorical(density.obs['radio'], categories = ['N', 'PNS', 'PS', 'S'])
        density = density[:,~density.var.index.isin(['Low Expr.', 'Other'])]

        for cond in metadata['CONDITIONS']:
            for pval_form in ['star', 'sci_notation']:
                if f'{cond}_color' in metadata:
                    palette = metadata[f'{cond}_color']
                else:
                    palette = None #'tab10'

                plot_grouped_key_mwu(
                    density,
                    condition_keys = [cond],
                    palette = palette,
                    pval_form = pval_form,
                    save_dir = f'figures/{p}_{celltype}_density_across_condition_filtered/'
                )


conditions = metadata['CONDITIONS'] + ['PID','description']
'''
Perform three comparisons at a patient level

1. Between tumor samples
2. Between normal and tumor samples combined. Here we take an approach saying that 
    normal samples and tumor samples are completely independent identity, and not from
    the same patient
3. Tumor - normal differentials per patient
'''
for p in ['PANEL_G', 'PANEL_H']:

    for celltype in metadata['CELLTYPES']:
        normal = sc.read(metadata[p]['AnnData'][f'patient_normal_{celltype}_name'])
        tumor = sc.read(metadata[p]['AnnData'][f'patient_tumor_{celltype}_name'])

        normal.X = normal.X / np.array(normal.obs['ROI_area'])[:,None]
        tumor.X = tumor.X / np.array(tumor.obs['ROI_area'])[:,None]

        n = normal.copy()
        t = tumor.copy()
        n.obs.index = n.obs.index + '_N'
        t.obs.index = t.obs.index + '_T'
        
        n2 = normal.copy()
        t2 = tumor.copy()
        intersection = n2.obs.index.intersection(t2.obs.index)
        n2 = n2[intersection]
        t2 = t2[intersection]
        sc.pp.log1p(n2)
        sc.pp.log1p(t2)
        t2.X = t2.X - n2.X

        comps = {
            'tumor': tumor,
            'combined': anndata.concat([n, t], axis = 0),
            'differential': t2,
        }

        os.makedirs('results/patient_csvs/', exist_ok = True)
        comps['combined'].write_csvs(f'results/patient_csvs/{p}_{celltype}_counts', skip_data = False)
        for c in comps:
            density = comps[c]
            for cond in density.obs.columns:
                density.obs[cond] = pd.Categorical(density.obs[cond])
            for pval_form in ['star', 'sci_notation']:
                if f'{cond}_color' in metadata:
                    palette = metadata[f'{cond}_color']
                else:
                    palette = None #'tab10'

                plot_grouped_key_mwu(
                    density,
                    condition_keys = metadata['CONDITIONS'],
                    palette = palette,
                    pval_form = pval_form,
                    save_dir = f'figures/patient_level_{p}_{celltype}_{c}_density_across_condition/'
                )


for p in ['PANEL_G', 'PANEL_H']:

    for celltype in metadata['CELLTYPES']:
        normal = sc.read(metadata[p]['AnnData'][f'patient_normal_{celltype}_name'])
        tumor = sc.read(metadata[p]['AnnData'][f'patient_tumor_{celltype}_name'])

        n = normal.copy()
        t = tumor.copy()
        n.obs.index = n.obs.index + '_N'
        t.obs.index = t.obs.index + '_T'
        
        n2 = normal.copy()
        t2 = tumor.copy()
        intersection = n2.obs.index.intersection(t2.obs.index)
        n2 = n2[intersection]
        t2 = t2[intersection]
        sc.pp.log1p(n2)
        sc.pp.log1p(t2)
        t2.X = t2.X - n2.X

        comps = {
            'tumor': tumor,
            'combined': anndata.concat([n, t], axis = 0),
            'differential': t2,
        }

        filt = [{'pathology': 'Normal', 'radio': 'N'}, {'pathology': 'AIS', 'radio': 'PNS'}, {'pathology': 'MIA', 'radio': 'PS'}, {'pathology': 'IAC', 'radio': 'S'}]


        for c in comps:
            density = comps[c]
            print(c)
            density = anndata.concat([
                density[(density.obs['pathology'] == s['pathology']) & (density.obs['radio'] == s['radio'])] for s in filt
            ])
            
            if 'N' in density.obs['radio']:
                density.obs['pathology'] = pd.Categorical(density.obs['pathology'], categories = ['Normal', 'AIS', 'MIA', 'IAC'])
                density.obs['radio'] = pd.Categorical(density.obs['radio'], categories = ['N', 'PNS', 'PS', 'S'])
            else:
                density.obs['pathology'] = pd.Categorical(density.obs['pathology'], categories = ['AIS', 'MIA', 'IAC'])
                density.obs['radio'] = pd.Categorical(density.obs['radio'], categories = ['PNS', 'PS', 'S'])
            
            for pval_form in ['star', 'sci_notation']:
                if f'{cond}_color' in metadata:
                    palette = metadata[f'{cond}_color'][:len(density.obs[cond].unique())]
                else:
                    palette = None #'tab10'

                plot_grouped_key_mwu(
                    density,
                    condition_keys = metadata['CONDITIONS'],
                    palette = palette,
                    pval_form = pval_form,
                    save_dir = f'figures/patient_level_{p}_{celltype}_{c}_density_across_condition/filtered/'
                )

        combined = anndata.concat([normal, tumor], axis = 0)
        combined = anndata.concat([
            combined[(combined.obs['pathology'] == s['pathology']) & (combined.obs['radio'] == s['radio'])] for s in filt
        ])
        sc.pp.log1p(combined)
        sc.pp.scale(combined)
        sc.pp.pca(combined)
        from sklearn.cluster import AgglomerativeClustering
        X_pca = combined.obsm['X_pca']
        cluster = AgglomerativeClustering(n_clusters=6, affinity='euclidean', linkage='ward')
        combined.obs['hclust_6'] = cluster.fit_predict(X_pca).astype(str)
        sc.pp.neighbors(combined, n_neighbors = 4)
        sc.tl.umap(combined)
        sc.pl.umap(combined, color = 'hclust_6', show = False)

        for cond in ['radio', 'pathology']:
            sc.pl.heatmap(combined, groupby = cond, var_names = combined.var.index, standard_scale = 'var', swap_axes = True, cmap = 'RdBu_r', dendrogram = True, save = f'Patient_{p}_{celltype}_heatmap_{cond}_filtered.pdf', show = False)
    

for p in ['PANEL_G', 'PANEL_H']:

    for celltype in metadata['CELLTYPES']:
        normal = sc.read(metadata[p]['AnnData'][f'patient_normal_{celltype}_name'])
        tumor = sc.read(metadata[p]['AnnData'][f'patient_tumor_{celltype}_name'])
        normal.X = normal.X / np.array(normal.obs['ROI_area'])[:,None]
        tumor.X = tumor.X / np.array(tumor.obs['ROI_area'])[:,None]
        normal.obs.index = normal.obs.index + '_N'
        tumor.obs.index = tumor.obs.index + '_T'
        
        combined = anndata.concat([normal, tumor], axis = 0)
        sc.pp.log1p(combined)
        sc.pp.scale(combined)
        sc.pp.pca(combined)
        from sklearn.cluster import AgglomerativeClustering
        X_pca = combined.obsm['X_pca']
        cluster = AgglomerativeClustering(n_clusters=6, affinity='euclidean', linkage='ward')
        combined.obs['hclust_6'] = cluster.fit_predict(X_pca).astype(str)
        sc.pp.neighbors(combined, n_neighbors = 4)
        sc.tl.umap(combined)
        sc.pl.umap(combined, color = 'hclust_6', show = False)

        for cond in ['radio', 'pathology']:
            sc.pl.heatmap(combined, groupby = cond, var_names = combined.var.index, standard_scale = 'var', swap_axes = True, cmap = 'RdBu_r', dendrogram = True, save = f'Patient_{p}_{celltype}_heatmap_{cond}.pdf', show = False)
        
        plt.close()





for p in ['PANEL_G', 'PANEL_H']:

    adata = adata_dict[p]
    adata.obs['PID'] = adata.obs['description'].str[:6]

    for celltype in metadata['CELLTYPES']: 
        normal = sc.read(metadata[p]['AnnData'][f'patient_normal_{celltype}_name'])
        tumor = sc.read(metadata[p]['AnnData'][f'patient_tumor_{celltype}_name'])
        normal.X = normal.X / np.array(normal.obs['ROI_area'])[:,None]
        tumor.X = tumor.X / np.array(tumor.obs['ROI_area'])[:,None]
        normal.obs.index = normal.obs.index + '_N'
        tumor.obs.index = tumor.obs.index
        
        combined = anndata.concat([normal, tumor], axis = 0)

        filt = [{'pathology': 'Normal', 'radio': 'N'}, {'pathology': 'AIS', 'radio': 'PNS'}, {'pathology': 'MIA', 'radio': 'PS'}, {'pathology': 'IAC', 'radio': 'S'}]

        a = grouped_obs_mean(adata[adata.obs['radio']!='N'], group_key = 'PID')
        ad = anndata.AnnData(a.loc[:,~a.columns.str.contains('EMPTY')], obs = tumor.obs)
        ad = anndata.concat([
            ad[(ad.obs['pathology'] == s['pathology']) & (ad.obs['radio'] == s['radio'])] for s in filt
        ])
        sc.pp.scale(ad, max_value = 2)
        ad.X[ad.X < -2] = -2
        for cond in ['radio', 'pathology']:
            sc.pl.heatmap(ad, groupby = cond, var_names = ad.var.index, standard_scale = 'obs', swap_axes = True, cmap = 'RdBu_r', dendrogram = True, save = f'Patient_{p}_{celltype}_heatmap_{cond}_mean_expression.pdf', show = False)
        


# Comparison between Normal - Tumor axis
# No evidence to support the claim
for p in ['PANEL_G', 'PANEL_H']:

    for celltype in metadata['CELLTYPES']:
        normal = sc.read(metadata[p]['AnnData'][f'patient_normal_{celltype}_name'])
        tumor = sc.read(metadata[p]['AnnData'][f'patient_tumor_{celltype}_name'])

        normal.X = normal.X / np.array(normal.obs['ROI_area'])[:,None]
        tumor.X = tumor.X / np.array(tumor.obs['ROI_area'])[:,None]

        n = normal.copy()
        t = tumor.copy()
        n.obs.index = n.obs.index + '_N'
        t.obs.index = t.obs.index + '_T'
        
        n2 = normal.copy()
        t2 = tumor.copy()
        intersection = n2.obs.index.intersection(t2.obs.index)
        n2 = n2[intersection]
        t2 = t2[intersection]

        n2.obs['radio'] = t2.obs['radio']
        n2.obs['pathology'] = t2.obs['pathology']

        filt = [{'pathology': 'Normal', 'radio': 'N'}, {'pathology': 'AIS', 'radio': 'PNS'}, {'pathology': 'MIA', 'radio': 'PS'}, {'pathology': 'IAC', 'radio': 'S'}]
        
        n2 = anndata.concat([
            n2[(n2.obs['pathology'] == s['pathology']) & (n2.obs['radio'] == s['radio'])] for s in filt
        ])
        n2.obs['pathology'] = pd.Categorical(n2.obs['pathology'], categories = ['AIS', 'MIA', 'IAC'])
        n2.obs['radio'] = pd.Categorical(n2.obs['radio'], categories = ['PNS', 'PS', 'S'])

        for cond in metadata['CONDITIONS']:
            for pval_form in ['star', 'sci_notation']:
                if f'{cond}_color' in metadata:
                    palette = metadata[f'{cond}_color']
                else:
                    palette = None #'tab10'

                plot_grouped_key_mwu(
                    n2,
                    condition_keys = [cond],
                    palette = palette,
                    pval_form = pval_form,
                    save_dir = f'figures/{p}_{celltype}_density_across_condition_normal_tumor_axis/'
                )


# patient feature selection
# measure mean intensity on patient level
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


conditions = metadata['CONDITIONS'] + ['PID', 'Race', 'Smoking', 'Packyrs', 'Age', 'Molecular']
for p in ['PANEL_G', 'PANEL_H']:

    adata = adata_dict[p]
    adata.obs['PID'] = adata.obs['description'].str[:6]

    n = adata[adata.obs['radio']=='N']
    t = adata[adata.obs['radio']!='N']

    n_ = grouped_obs_mean(n, group_key = 'PID')
    n_obs = n.obs[conditions].drop_duplicates().sort_values('PID').set_index('PID')
    n_obs.index.name = None
    n_data = anndata.AnnData(X = n_, obs = n_obs)
    n_data.obs['group'] = 'N'
    n_data = n_data[:,~n_data.var.index.str.contains('EMPTY')]
    n_data.obs.index = n_data.obs.index + '_N'
    
    t_ = grouped_obs_mean(t, group_key = 'PID')
    t_obs = t.obs[conditions].drop_duplicates().sort_values('PID').set_index('PID')
    t_obs.index.name = None
    t_data = anndata.AnnData(X = t_, obs = t_obs)
    t_data.obs['group'] = 'T'
    t_data = t_data[:,~t_data.var.index.str.contains('EMPTY')]
    t_data.obs.index = t_data.obs.index + '_T'

    n_data.write(metadata[p]['AnnData']['patient_normal_expression_name'])
    t_data.write(metadata[p]['AnnData']['patient_tumor_expression_name'])



normal.var['group'] = 'N'
tumor.var['group'] = 'T'
adata = anndata.concat([normal, tumor], axis = 1, join = 'outer')
adata.obs['pathology'] = tumor.obs['pathology']
adata.obs['radio'] = tumor.obs['radio']
adata.obs['ROI_area_T'] = tumor.obs['ROI_area']
adata.obs['ROI_area_N'] = normal.obs['ROI_area']


adata.write('results/Patient_cellcount_all.h5ad')



# # beware of NaNs
# sc.pp.scale(adata)
# sc.pp.pca(adata)
# from sklearn.cluster import AgglomerativeClustering
# cluster = AgglomerativeClustering(n_clusters=5, affinity='euclidean', linkage='ward')
# adata.obs['hclust_5'] = cluster.fit_predict(X_pca).astype(str)