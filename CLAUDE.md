# ggo-imc — Claude Code Configuration

## Sync Before Work

1. @Mission.md — current todo list and phase status
2. @journal_summary.md — recent decisions

For repo-wide rules (container, conventions, testing): see repo root CLAUDE.md.

---

## Project Context

**Goal:** IMC analysis of GGO lung samples — cell type characterization, differential abundance, spatial microenvironment, patient-level analysis.

**Platform:** Imaging Mass Cytometry (IMC)

**Note on R:** The `patient` Snakemake rule runs `scripts/asd.R` directly (no container) — the sc_tools image does not include R. Run R rules in the `ggo_imc` conda env with R and dependencies installed separately.

---

## Running This Project

```bash
# From repo root:
./scripts/run_container.sh projects/imc/ggo-imc python scripts/<script>.py
snakemake -d projects/imc/ggo-imc -s projects/imc/ggo-imc/Snakefile <target>

# Local conda (no container):
conda activate ggo_imc
SC_TOOLS_RUNTIME=none snakemake -d projects/imc/ggo-imc -s projects/imc/ggo-imc/Snakefile <target>

# Tests:
pytest projects/imc/ggo-imc/tests/ -v
```

**Snakemake targets (by figure):** `figure1`, `figure2`, `figure3`, `figure4`, `figure5`, `all`

---

## Key Files

| Path | Description |
|------|-------------|
| `scripts/celltype_heatmap_info.py` | Figure 1: cell type heatmap |
| `scripts/celltype_differential_abundance.py` | Figure 2/3: immune/stromal density |
| `scripts/t_cell_analysis.py` | Figure 2: T cell analysis |
| `scripts/myeloid_analysis.py` | Figure 2: myeloid analysis |
| `scripts/ue_analysis.py` | Figure 4: microenvironment |
| `scripts/patient_group.py` | Figure 5: patient risk groups |
| `scripts/asd.R` | Figure 5: patient-level R analysis (runs outside container) |
