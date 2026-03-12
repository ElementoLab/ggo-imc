import imc_analysis as imc
import scanpy as sc
import sys
import pandas as pd
# import matplotlib
# matplotlib.use('Agg')

metadata = imc.utils.parse_yaml('metadata/ggo_config.yml')

print(sys.argv)

# reading in data
adata_dict = dict()
for p in ['PANEL_G', 'PANEL_H']:
    print(f"Reading {metadata[p]['AnnData']['phenotyped_umap_name']}...")
    adata_dict[p] = sc.read(
        metadata[p]['AnnData']['phenotyped_umap_name'],
        backup_url = metadata[p]['AnnData']['backup_url'])
    adata_dict[p].obs['Smoking Status'] = adata_dict[p].obs['Smoking Status'].astype(str).str.replace(' smokes', '')
    adata_dict[p].obs['Smoking Status'] = pd.Categorical(adata_dict[p].obs['Smoking Status'],
        categories = ['Never', 'Former', 'Current'], ordered = True)
    adata_dict[p].uns['Smoking Status_colors'] = ['#ffffff', '#9D9086', '#725b43']


cond ='Smoking Status'
# differential cell type densities
for p in ['PANEL_G', 'PANEL_H']:
    adata = adata_dict[p][adata_dict[p].obs['pathology'] != 'N'].copy()
    
    # computing celltype density
    density = imc.tl.celltype_density(
        adata,
        celltype = 'celltype_broad',
        condition_keys = ['Smoking Status'])

    # statistical testing
    imc.tl.grouped_mwu_test(
        density,
        condition_keys = [cond]
    )

    for pval_form in ['star', 'sci_not']:
        imc.pl.plot_mwu(
            density,
            kind = 'box-line',
            save_dir=f'figures/densities_{cond}_tumor/{p}/',
            pval_form=pval_form
        )


    adata = adata_dict[p][adata_dict[p].obs['pathology'] == 'N'].copy()
    
    # computing celltype density
    density = imc.tl.celltype_density(
        adata,
        celltype = 'celltype_broad',
        condition_keys = ['Smoking Status'])

    # statistical testing
    imc.tl.grouped_mwu_test(
        density,
        condition_keys = [cond]
    )

    for pval_form in ['star', 'sci_not']:
        imc.pl.plot_mwu(
            density,
            kind = 'box-line',
            save_dir=f'figures/densities_{cond}_normal/{p}/',
            pval_form=pval_form
        )