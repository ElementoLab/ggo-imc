# Phase 1: Analytical Robustness - Context

**Gathered:** 2026-03-20
**Status:** Ready for planning

<domain>
## Phase Boundary

All new computational analyses are complete, validated, and ready to feed into figures and text. This phase produces: clustering validation outputs (silhouette + alternative clustering), G4 diagnostic performance metrics, Yoffe cohort stratification, and clinical covariate breakdowns by G1–G4. No figure display work — outputs are validated results and supplementary-ready plots.

</domain>

<decisions>
## Implementation Decisions

### Clustering Validation (CLUS-01, CLUS-02)
- Silhouette analysis over k=2 to 6 (tight range centered on k=4)
- Output: mean silhouette width line plot (k on x-axis, mean silhouette width on y-axis) — one clean supplementary figure
- Alternative method: k-means only (k=4 to mirror hierarchical); no Leiden
- Concordance reported as ARI + % agreement — both numbers, fits in one sentence of reviewer response text

### G4 Sensitivity/Specificity (CLUS-03)
- Positive class: radiologically pure non-solid (PNS) patients — main mode of detection, least invasive
- Operationalization: G4 membership as predictor of PNS histology; sensitivity = PNS cases in G4 / all PNS cases; specificity = non-PNS cases not in G4 / all non-PNS cases
- Report with Wilson or Clopper-Pearson 95% CIs (note: n=8 PNS cases → CIs will be wide)
- Presentation: inline text in Results, one sentence; no new figure panel

### Yoffe Cohort Stratification (ANA-01)
- Method: centroid-based nearest neighbor — compute G1–G4 centroids in process-score space from the 122 IMC patients, assign each Yoffe patient to nearest centroid
- Feature space: IMC process scores (EMT, B&T, Fibrosis, Angiogenesis, Epithelial, Mesenchymal, Pan-Immune) — same scores computed in patient_group.py, generatable from Yoffe spatial transcriptomics using the same cell-type gene lists
- Success criterion: G4-enriched in high-risk histology (IAC/MIA) in Yoffe vs AAH/AIS — tested with chi-squared or Fisher's exact

### Clinical Covariate Breakdowns (ANA-02, ANA-03, ANA-04, ANA-05)
- Covariates: driver mutations (EGFR, KRAS, RBM10, TP53), smoking status, IASLC grading, PNS/IAC case counts
- Statistical tests: chi-squared or Fisher's exact per covariate (Fisher's when any cell count <5)
- Output format: stacked bar plots per covariate in supplementary — one plot per covariate showing % per group with p-value annotation
- PNS/IAC case counts (ANA-05): compute in Phase 1; Phase 2 folds into Fig 5d/5e display

### Adjusted P-values and Inter-panel Consistency (ANA-06, ANA-07)
- ANA-06: Confirm and report adjusted p-values for Fig 2l–m M1 macrophage N vs AIS comparisons
- ANA-07: Document inter-panel consistency between tumor and immune IMC panels in analysis record (feeds methods text in Phase 3)

### Claude's Discretion
- Exact silhouette distance metric (euclidean vs cosine) — choose based on data dimensionality
- k-means random seed and n_init — standard settings
- Supplementary figure layout details
- Whether to write a new script or extend patient_group.py for clustering validation

</decisions>

<canonical_refs>
## Canonical References

**Downstream agents MUST read these before planning or implementing.**

### Phase requirements and success criteria
- `.planning/ROADMAP.md` — Phase 1 goal, success criteria (5 items), and requirement traceability
- `.planning/REQUIREMENTS.md` — CLUS-01, CLUS-02, CLUS-03, ANA-01 through ANA-07 requirement specs

### Existing pipeline code
- `scripts/patient_group.py` — Fig 5 pipeline; loads patient_celltype_broad_clustered AnnData, computes process scores via sc.tl.score_genes, has PCA/neighbors/UMAP infrastructure; new clustering analyses should extend or parallel this script
- `scripts/utils.py` — shared load_config, save_figure, ensure_dir; all new scripts must use these
- `scripts/exploratory/mutation.py` — existing mutation analysis (reference for clinical covariate work)
- `scripts/exploratory/clinical_features.py` — existing clinical feature analysis (reference for covariate breakdowns)

### Runtime constraint
- `scripts/asd.R` — Fig 5 patient analysis runs outside container in `ggo_imc` conda env; any R-based analyses must account for this

</canonical_refs>

<code_context>
## Existing Code Insights

### Reusable Assets
- `scripts/patient_group.py`: Already loads the correct AnnData (`patient_celltype_broad_clustered`), computes all 7 process scores via `sc.tl.score_genes`, and has scanpy PCA/neighbors/UMAP in place. Silhouette and k-means analyses should load the same AnnData and use the same scaled/PCA representation.
- `scripts/utils.py`: `load_config()`, `save_figure()`, `ensure_dir()` — use in all new scripts
- `scripts/exploratory/mutation.py` and `clinical_features.py`: Reference these for how metadata is loaded and joined; avoid duplicating that logic

### Established Patterns
- Config-driven: all paths come from `metadata/ggo_config.yml` via `load_config()`
- Figures saved via `save_figure()` to `figures/<subfolder>/`
- AnnData `.obs` contains patient-level metadata including `Group`, `radio`, `pathology` columns

### Integration Points
- Phase 1 outputs feed Phase 2 figure scripts directly: clinical covariate plots → supplementary; PNS/IAC counts → Fig 5d/5e; silhouette plot → supplementary
- Yoffe cohort AnnData already integrated (per PROJECT.md) — confirm path in ggo_config.yml before planning

</code_context>

<specifics>
## Specific Ideas

- G4 molecular profile: low immune + high stroma (fibroblast/myofibroblast) — this is the biological rationale for the group, not just a label; sens/spec framing should make this clear in the reviewer response
- PNS patients are the key clinical anchor for G4: radiological PNS = least-invasive detection mode; G4 capturing PNS patients = clinical relevance of the molecular signature

</specifics>

<deferred>
## Deferred Ideas

None — discussion stayed within phase scope.

</deferred>

---

*Phase: 01-analytical-robustness*
*Context gathered: 2026-03-20*
