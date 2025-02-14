import imc_analysis as imc
import scanpy as sc

import os
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

'''
Figure 5 Patient PCA
'''
# color by Patient group
a = sc.read(
    metadata['patient_celltype_broad_clustered'],
    backup_url = metadata['patient_group_url']
)

# concatenate patient group data
tmp = density_all.obs.merge(a.obs[['Group']], left_on = 'GGO ID', right_index = True)
tmp['Group'] = tmp['Group'].astype(str)
tmp.loc[tmp['pathology']=='N','Group'] = 'N'
tmp['Group'] = pd.Categorical(tmp['Group'], categories = ['N', 'Group1', 'Group2', 'Group3', 'Group4'], ordered = True)
da = density_all[tmp.index]
da.obs = tmp
da.uns['Group_colors'] = ['#2AC674', '#8DAC91', '#F9C8AB', '#DB4962', '#252F32']


# Plot PCA
fig, ax = plt.subplots(1,1,figsize = (6,5), dpi = 300)
sc.pl.pca(
    da,
    color = 'Group',
    size = 200,
    add_outline = True,
    outline_width = (0.15, 0.),
    annotate_var_explained = True,
    title = '',
    legend_loc = None,
    ax = ax,
    show = False
) 

# display mean and covariance matrix as ellipse on pca plot
for i,group in enumerate(da.obs['Group'].cat.categories):
    da_group = da[da.obs['Group'] == group]
    mean = da_group.obsm['X_pca'][:,:2].mean(axis = 0)
    cov = np.cov(da_group.obsm['X_pca'][:,:2].T)

    from matplotlib.patches import Ellipse

    eigenvalues, eigenvectors = np.linalg.eig(cov)
    angle = np.degrees(np.arctan2(*eigenvectors[:, 0][::-1]))
    width, height = 2 * np.sqrt(eigenvalues)
    ellipse = Ellipse(xy=mean, width=width, height=height, angle=angle, facecolor=da.uns['Group_colors'][i], edgecolor = None, alpha = 0.3)
    plt.gca().add_patch(ellipse)

os.makedirs('figure/figure5', exist_ok = True)
sns.despine()
plt.tight_layout()
plt.savefig('figures/figure5/roi_pca_group.pdf', bbox_inches = 'tight')
plt.close()

