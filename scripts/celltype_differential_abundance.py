import imc_analysis as imc
import scanpy as sc
import sys
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

# differential cell type densities
for p in ['PANEL_G', 'PANEL_H']:
    # computing celltype density
    density = imc.tl.celltype_density(
        adata_dict[p],
        celltype = 'celltype_broad',
        condition_keys = metadata['CONDITIONS'])

    # per condition
    for cond in ['pathology', 'radio']:
        # statistical testing
        imc.tl.grouped_mwu_test(
            density,
            condition_keys = [cond]
        )

        densities = dict()
        # if lineage is provided as sys.arg, subset cell lineage
        lineages = metadata['CELL_LINEAGES']
        x = set(sys.argv).intersection(set(lineages.keys()))
        
        if len(x):
            print(f'Subsetting samples for cell lineage: "{x}"')

            cell_lineages = dict()
            for x in lineages.keys():
                if x in sys.argv:
                    cell_lineages[x] = lineages[x]
            for lineage in cell_lineages:
                densities[lineage] = density[:,density.var.index.intersection(cell_lineages[lineage])]
        else:
            densities['all cells'] = density

        # produce figures
        for lineage in densities:

            if lineage in ['lymphoid', 'myeloid']:
                figure_number = 'figure2'
            elif lineage in ['stromal', 'epithelial']:
                figure_number = 'figure3'
            else:
                figure_number = ''

            for pval_form in ['star', 'sci_notation']:
                imc.pl.plot_mwu(
                    densities[lineage],
                    kind = 'box-line',
                    save_dir=f'figures/{figure_number}/densities_{lineage}/{p}/',
                    pval_form=pval_form
                )
