
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

import yaml
import matplotlib
sc.settings.set_figure_params(dpi=200, dpi_save=300, fontsize=12)
matplotlib.rcParams["pdf.fonttype"] = 42
matplotlib.rcParams["ps.fonttype"] = 42
matplotlib.rcParams["axes.grid"] = False


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

#matplotlib.use('Agg')
metadata_filename = 'metadata/liver_config.yml'

with open(metadata_filename, "r") as stream:
    try:
        metadata = yaml.safe_load(stream)
    except yaml.YAMLError as exc:
        print(exc)


for celltype in metadata['CELLTYPES']:

    patient = sc.read(metadata[f'patient_{celltype}_name'])

    patient_density = patient.copy()
    patient_density.X = patient_density.X / np.array(patient.obs['ROI_area'])[:,None]
    for cond in patient_density.obs:
        patient_density.obs[cond] = pd.Categorical(patient_density.obs[cond])

    for pval_form in ['star', 'sci_notation']:
        if f'{cond}_color' in metadata:
            palette = metadata[f'{cond}_color']
        else:
            palette = None #'tab10'

        plot_grouped_key_mwu(
            patient_density,
            condition_keys = metadata['CONDITIONS'],
            palette = palette,
            pval_form = pval_form,
            save_dir = f'figures/patient_{celltype}_density_across_condition/'
        )
