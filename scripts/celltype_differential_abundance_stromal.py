import imc_analysis as imc
import scanpy as sc

# import matplotlib
# matplotlib.use('Agg')

metadata = imc.utils.parse_yaml('metadata/ggo_config.yml')

# reading in data
adata_dict = dict()
for p in ['PANEL_G', 'PANEL_H']:
    print(f"Reading {metadata[p]['AnnData']['phenotyped_umap_name']}...")
    adata_dict[p] = sc.read(
        metadata[p]['AnnData']['phenotyped_umap_name'],
        backup_url = metadata[p]['AnnData']['backup_url'])

# differential cell type densities
celltype = 'celltype_broad'
for p in ['PANEL_G', 'PANEL_H']:
    # computing celltype density
    density = imc.tl.celltype_density(
        adata_dict[p],
        celltype = celltype,
        condition_keys = metadata['CONDITIONS'])

    density.write(metadata[p]['AnnData'][f'roi_{celltype}_name'])

    # per condition
    for cond in ['pathology', 'radio']:
        # statistical testing
        imc.tl.grouped_mwu_test(
            density,
            condition_keys = [cond]
        )

        # produce figures
        for pval_form in ['star', 'sci_notation']:
            imc.pl.plot_mwu(
                density,
                kind = 'box-line',
                save_dir=f'figures/figure2/differential_{celltype}_density/{p}/',
                pval_form=pval_form
            )
