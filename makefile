.PHONY: download celltype diff_abundance t_cell myeloid epithelial microenvironment pca patient # spatial similarity utag

download:
	python scripts/download_yaml.py

celltype:
	python scripts/celltype_heatmap_info.py

diff_abundance:
	python scripts/celltype_differential_abundance.py

t_cell:
	python scripts/t_cell_analysis.py

myeloid:
	python scripts/myeloid_analysis.py
	
epithelial:
	python scripts/epithelial_characterization.py

pca:
	python scripts/roi_pca_plot.py

pca_group:
	python scripts/roi_pca_plot_group.py

# spatial:
# 	python scripts/spatial_plot.py

# interaction:
# 	python scripts/cell_interaction.py

# similarity:
# 	python scripts/sample_similarity.py

microenvironment:
	python scripts/ue_analysis.py

patient:
	python scripts/create_patient_density_matrix.py
	r scripts/asd.R
	python scripts/patient_group.py

clean:
	$(RM) __pycache__
	$(RM) $(VENV)

run: setup download
figure1: celltype pca
figure2: diff_abundance
figure3: microenvironment
figure4: similarity
figure5: patient pca_group