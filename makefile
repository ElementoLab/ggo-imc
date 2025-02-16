.PHONY: download celltype diff_abundance t_cell myeloid epithelial microenvironment pca patient # spatial similarity utag

download:
	python scripts/download_yaml.py

# figure 1
celltype:
	python scripts/celltype_heatmap_info.py

pca:
	python scripts/roi_pca_plot.py

# figure 2
densities_immune:
	python scripts/celltype_differential_abundance.py lymphoid myeloid

t_cell:
	python scripts/t_cell_analysis.py

myeloid:
	python scripts/myeloid_analysis.py

# figure 3
densities_stromal_and_epithelial:
	python scripts/celltype_differential_abundance.py stromal epithelial

epithelial:
	python scripts/epithelial_characterization.py

# figure 4
microenvironment:
	python scripts/ue_analysis.py

# figure 5
pca_group:
	python scripts/roi_pca_plot_group.py

patient:
	r scripts/asd.R

patient_risk:
	python scripts/patient_group.py

# spatial:
# 	python scripts/spatial_plot.py

# interaction:
# 	python scripts/cell_interaction.py

# similarity:
# 	python scripts/sample_similarity.py

clean:
	$(RM) __pycache__

run: setup download
figure1: celltype pca
figure2: densities_immune t_cell myeloid
figure3: densities_stromal_and_epithelial epithelial
figure4: microenvironment
figure5: patient_risk pca_group
#figure5: patient patient_risk pca_group