import seaborn as sns
import matplotlib.pyplot as plt
import scanpy as sc
from tqdm import tqdm

pg = sc.read('results/panel_g.phenotyped.h5ad')
pg_epi = pg[pg.obs['celltype_broad'].isin(['Epi. Prol.', 'Epi.-like','Mesen.-like'])]
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
    x = 'PanCytokeratin(Pt198)', y = 'Vimentin(Nd143)', fill = True, ax = ax, cmap = 'YlOrRd')
    ax.set_title(rad)
    #ax.legend().remove()

plt.tight_layout()
plt.savefig(f'figures/EMT proportion_2.pdf', bbox_inches = 'tight')
plt.close()



