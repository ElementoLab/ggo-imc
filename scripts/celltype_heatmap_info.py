import imc_analysis as imc
import os
import scanpy as sc
import matplotlib.pyplot as plt

import matplotlib
matplotlib.use('Agg')

url = 'https://wcm.box.com/shared/static/mdntp2xf9tjobxeidkw93mg8jysb7nh9.yml'
metadata_path = 'metadata/ggo_config.yml'

# if metadata not found, download
imc.utils.download(url, metadata_path)

metadata = imc.utils.parse_yaml(metadata_path)

adata_dict = dict()
for p in ['PANEL_G', 'PANEL_H']:
    print(f"Reading {metadata[p]['AnnData']['phenotyped_umap_name']}...")
    adata_dict[p] = sc.read(
        metadata[p]['AnnData']['phenotyped_umap_name'],
        backup_url = metadata[p]['AnnData']['backup_url'])

# for celltype in celltypes:
#     adata_dict[p] = adata_dict[p][adata_dict[p].obs[celltype] != 'Low Expr.']

os.makedirs('figures/figure1/', exist_ok = True)
print('Plot celltype heatmap')
for p in adata_dict:
    imc.pl.celltype_heatmap(
        adata_dict[p],
        cluster_ids = ['celltype_broad'],
        var_names = metadata[p]['var_celltype_groups'],
        out_dir = f'figures/figure1/celltype/{p}/'
    )


print('Plot UMAP')
for p in adata_dict:
    adata = sc.pp.subsample(
        adata_dict[p],
        copy = True,
        n_obs = 100000,
        random_state = 0,
    )
    
    celltype = 'celltype_broad'

    fig, ax = plt.subplots(1,1,figsize = (6,5), dpi = 300)
    sc.pl.umap(adata, color = celltype, frameon = False, show = False)

    os.makedirs(f'figures/figure1/celltype/{p}/', exist_ok = True)
    plt.tight_layout()
    plt.savefig(f'figures/figure1/celltype/{p}/umap.pdf')
    plt.close()
