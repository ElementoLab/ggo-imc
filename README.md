# ggo-imc

This repository contains the code to reproduce figures from manuscript: "Simultaneous immunomodulation and epithelial-to-mesenchymal transition drives lung adenocarcinoma progression."

## Data
The processed data and raw data (for further exploration) can be found at Zenodo (view of the data might be restricted until publication):

- [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.14822106.svg)](https://doi.org/10.5281/zenodo.14822106) - Processed Data
- [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.14822106.svg)](https://doi.org/10.5281/zenodo.14822106) - Raw Data

## Environment Setup
```
# setting up conda environment
conda create -n ggo-imc python==3.9 -y
conda activate ggo-imc
pip install -r requirements.txt

# downloading data
make download
```

## Analysis Pipeline
For streamlined analysis to produce figures use the command
```
make run
```

Alternatively for individual plots:
```
make figure1
make figure2
make figure3
make figure4
make figure5
```

<!-- ```
# plot celltype heatmaps
make celltype

# to plot differential abundance
make diff_abundance

# plot interactions
make interaction

make t_cell

make myelod

make epithelial

make microenvironment

``` -->
