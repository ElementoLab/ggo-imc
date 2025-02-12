import imc_analysis as imc
import scanpy as sc

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
sc.settings.set_figure_params(dpi=200, dpi_save=300, fontsize=12)
import matplotlib

# load data
metadata = imc.utils.parse_yaml('metadata/ggo_config.yml')

density_all = sc.read(metadata['PCA_ANNDATA_NAME'],
    backup_url = metadata['PCA_URL'])

# Perform PCA
# Flip x axis such that normal is on the negative side
sc.pp.scale(density_all)
sc.pp.pca(density_all)
density_all.obsm['X_pca'][:,0] = -density_all.obsm['X_pca'][:,0]

# Plot PCA
fig, ax = plt.subplots(1,1,figsize = (12,8), dpi = 300)
sc.pl.pca(
    density_all,
    color = 'pathology',
    size = 200,
    add_outline = True,
    outline_width = (0.15, 0.),
    annotate_var_explained = True,
    title = '',
    ax = ax,
    show = False
) 
sns.despine()
plt.tight_layout()
plt.savefig('figures/figure1/roi_pca_all_combined.pdf', bbox_inches = 'tight')
plt.close()
