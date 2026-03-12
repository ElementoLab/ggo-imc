# ggo-imc IMC Snakemake pipeline
# Run from repo root: snakemake -d projects/imc/ggo-imc -s projects/imc/ggo-imc/Snakefile [target]
# Or from project dir: snakemake -d . -s Snakefile [target]
# Runtime: auto-detected via run_container.sh (Docker on macOS, Apptainer on Linux).
# Local conda: conda activate ggo_imc && SC_TOOLS_RUNTIME=none snakemake -d . -s Snakefile <target>
# Note: the patient rule runs R (scripts/asd.R) directly — sc_tools image does not include R.

configfile: "config.yaml"

ROOT = config["repo_root"]
PROJECT = config["project_rel"]

def run_container(script, args=""):
    cmd = f'cd {ROOT} && ./scripts/run_container.sh {PROJECT} python {script}'
    return cmd + (f' {args}' if args else '')

# ---- Setup / download ----
rule download:
    output: touch("download.done")
    input: "scripts/download_yaml.py"
    shell: run_container("scripts/download_yaml.py") + " && touch download.done"

# ---- Figure 1 ----
rule celltype:
    output: touch("celltype.done")
    input: "scripts/celltype_heatmap_info.py"
    shell: run_container("scripts/celltype_heatmap_info.py") + " && touch celltype.done"

rule pca:
    output: touch("pca.done")
    input: "scripts/roi_pca_plot.py"
    shell: run_container("scripts/roi_pca_plot.py") + " && touch pca.done"

# ---- Figure 2 ----
rule densities_immune:
    output: touch("densities_immune.done")
    input: "scripts/celltype_differential_abundance.py"
    shell: run_container("scripts/celltype_differential_abundance.py", "lymphoid myeloid") + " && touch densities_immune.done"

rule t_cell:
    output: touch("t_cell.done")
    input: "scripts/t_cell_analysis.py"
    shell: run_container("scripts/t_cell_analysis.py") + " && touch t_cell.done"

rule myeloid:
    output: touch("myeloid.done")
    input: "scripts/myeloid_analysis.py"
    shell: run_container("scripts/myeloid_analysis.py") + " && touch myeloid.done"

# ---- Figure 3 ----
rule densities_stromal_and_epithelial:
    output: touch("densities_stromal_epithelial.done")
    input: "scripts/celltype_differential_abundance.py"
    shell: run_container("scripts/celltype_differential_abundance.py", "stromal epithelial") + " && touch densities_stromal_epithelial.done"

rule epithelial:
    output: touch("epithelial.done")
    input: "scripts/epithelial_characterization.py"
    shell: run_container("scripts/epithelial_characterization.py") + " && touch epithelial.done"

# ---- Figure 4 ----
rule microenvironment:
    output: touch("microenvironment.done")
    input: "scripts/ue_analysis.py"
    shell: run_container("scripts/ue_analysis.py") + " && touch microenvironment.done"

# ---- Figure 5 ----
rule pca_group:
    output: touch("pca_group.done")
    input: "scripts/roi_pca_plot_group.py"
    shell: run_container("scripts/roi_pca_plot_group.py") + " && touch pca_group.done"

rule patient:
    # NOTE: scripts/asd.R is not yet committed to this repo.
    # Run in ggo_imc conda env (R + dependencies); sc_tools image does not include R.
    output: touch("patient.done")
    input: "scripts/asd.R"
    shell: "Rscript scripts/asd.R && touch patient.done"

rule patient_risk:
    output: touch("patient_risk.done")
    input: "scripts/patient_group.py"
    shell: run_container("scripts/patient_group.py") + " && touch patient_risk.done"

# ---- Aggregates ----
rule figure1:
    input: "celltype.done", "pca.done"

rule figure2:
    input: "densities_immune.done", "t_cell.done", "myeloid.done"

rule figure3:
    input: "densities_stromal_epithelial.done", "epithelial.done"

rule figure4:
    input: "microenvironment.done"

rule figure5:
    input: "patient_risk.done", "pca_group.done"

rule run:
    input: "download.done"

rule all:
    input:
        "celltype.done", "pca.done",
        "densities_immune.done", "t_cell.done", "myeloid.done",
        "densities_stromal_epithelial.done", "epithelial.done",
        "microenvironment.done",
        "patient_risk.done", "pca_group.done"

rule clean:
    shell: "rm -f *.done"
