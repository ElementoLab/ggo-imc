.PHONY: download spatial celltype t_cell similarity utag diff_abundance pca patient

download:
	echo not implemented yet
	# to be comming

celltype:
	python scripts/celltype_heatmap_info.py

diff_abundance:
	python scripts/celltype_differential_abundance.py

t_cell:
	python scripts/t_cell_analysis.py

pca:
	python scripts/roi_pca_plot.py

spatial:
	python scripts/spatial_plot.py

interaction:
	python scripts/cell_interaction.py

similarity:
	python scripts/sample_similarity.py

utag:
	python scripts/utag_ue.py	

patient:
	python scripts/create_patient_density_matrix.py
	r scripts/asd.R
	python scripts/patient_group.py

clean:
	$(RM) __pycache__
	$(RM) $(VENV)

run: setup download
figure1: celltype
figure2: diff_abundance
figure3: diff_abundance
figure4: similarity
figure5: patient