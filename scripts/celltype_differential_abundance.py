"""Figures 2/3: differential cell type density by pathology and radiology.

Accepts optional positional arguments (cell lineage names from CELL_LINEAGES
in the config) to subset the output figures to specific lineages.
"""

import sys

import imc_analysis as imc
import scanpy as sc

from utils import load_config, load_panels


if __name__ == '__main__':
    metadata = load_config()
    print(sys.argv)

    adata_dict = load_panels(metadata)
    for p in ['PANEL_G', 'PANEL_H']:
        adata_dict[p].uns['pathology_colors'] = metadata['pathology_color']
        adata_dict[p].uns['radio_colors'] = metadata['Radiology_color']

    # differential cell type densities
    for p in ['PANEL_G', 'PANEL_H']:
        density = imc.tl.celltype_density(
            adata_dict[p],
            celltype='celltype_broad',
            condition_keys=metadata['CONDITIONS'],
        )

        for cond in ['pathology', 'radio']:
            imc.tl.grouped_mwu_test(density, condition_keys=[cond])

            densities: dict = {}
            # CELL_LINEAGES is optional in ggo_config.yml.
            # Default to an empty dict so the else-branch runs all cells.
            lineages = metadata.get('CELL_LINEAGES', {})
            x = set(sys.argv).intersection(set(lineages.keys()))

            if len(x):
                print(f'Subsetting samples for cell lineage: "{x}"')
                cell_lineages = {k: lineages[k] for k in lineages if k in sys.argv}
                for lineage, members in cell_lineages.items():
                    densities[lineage] = density[:, density.var.index.intersection(members)]
            else:
                densities['all cells'] = density

            for lineage, dens in densities.items():
                if lineage in ['lymphoid', 'myeloid']:
                    figure_number = 'figure2'
                elif lineage in ['stromal', 'epithelial']:
                    figure_number = 'figure3'
                else:
                    figure_number = ''

                for pval_form in ['star', 'sci_notation']:
                    imc.pl.plot_mwu(
                        dens,
                        kind='box-line',
                        save_dir=f'figures/{figure_number}/densities_{lineage}/{p}/',
                        pval_form=pval_form,
                    )
