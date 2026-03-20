# Phase 1: Analytical Robustness - Research

**Researched:** 2026-03-20
**Domain:** Computational biology / IMC single-cell analysis / statistical validation
**Confidence:** HIGH — all findings verified directly from codebase inspection and live data

---

<user_constraints>
## User Constraints (from CONTEXT.md)

### Locked Decisions
- Silhouette analysis over k=2 to 6 (tight range centered on k=4)
- Output: mean silhouette width line plot (k on x-axis, mean silhouette width on y-axis) — one clean supplementary figure
- Alternative method: k-means only (k=4 to mirror hierarchical); no Leiden
- Concordance reported as ARI + % agreement — both numbers, fits in one sentence of reviewer response text
- G4 sens/spec: Positive class is radiologically pure non-solid (PNS) patients
- Operationalization: G4 membership as predictor of PNS histology; sensitivity = PNS cases in G4 / all PNS cases; specificity = non-PNS cases not in G4 / all non-PNS cases
- Report with Wilson or Clopper-Pearson 95% CIs (n=8 PNS-IAC cases → CIs will be wide)
- Yoffe stratification: centroid-based nearest neighbor in process-score space (7 scores from patient_group.py)
- Yoffe success criterion: G4-enriched in high-risk histology (IAC/MIA) vs AAH/AIS — tested with chi-squared or Fisher's exact
- Clinical covariate output: stacked bar plots per covariate in supplementary
- Statistical tests: chi-squared or Fisher's exact (Fisher's when any cell count <5)
- ANA-06: Confirm adjusted p-values for Fig 2l–m M1 macrophage N vs AIS comparisons
- ANA-07: Document inter-panel consistency between tumor and immune IMC panels in analysis record

### Claude's Discretion
- Exact silhouette distance metric (euclidean vs cosine) — choose based on data dimensionality
- k-means random seed and n_init — standard settings
- Supplementary figure layout details
- Whether to write a new script or extend patient_group.py for clustering validation

### Deferred Ideas (OUT OF SCOPE)
None — discussion stayed within phase scope.
</user_constraints>

---

<phase_requirements>
## Phase Requirements

| ID | Description | Research Support |
|----|-------------|-----------------|
| CLUS-01 | Silhouette width analysis confirms k=4 as optimal patient grouping | AnnData `results/patient_celltype_broad_density_group.h5ad` has 122 patients × 14 cell types (X = densities); sklearn silhouette_score usable directly on scaled/PCA representation from patient_group.py |
| CLUS-02 | Alternative clustering method (k-means) applied and compared to hierarchical result | sklearn KMeans; Group labels already in adata.obs as 'Group1'–'Group4'; ARI via sklearn.metrics.adjusted_rand_score; % agreement via direct label comparison after remapping |
| CLUS-03 | Sensitivity/specificity of G4 for identifying PNS/IAC cases | adata.obs has 'radio' (PNS/PS/S) and 'Group'; n=61 PNS total; 8 PNS patients with IAC pathology; Wilson CI via statsmodels.stats.proportion or scipy.stats |
| ANA-01 | Yoffe et al. cohort stratified into G1–G4 | Yoffe AnnData path must be confirmed in ggo_config.yml (not currently present — key gap); centroid computation from process scores is straightforward with sklearn NearestCentroid or manual scipy.spatial.distance.cdist |
| ANA-02 | Driver mutation distribution (EGFR, KRAS, RBM10, TP53) across G1–G4 | adata.obs already has EGFR, KRAS, RBM10, TP53 as 'WT'/'MUT' categories (verified); chi-squared or Fisher's exact via scipy.stats |
| ANA-03 | Smoking status distribution across G1–G4 | adata.obs has 'Smoking Status'; categories: Never, Former smokes, Current smokes (verified in clinical_features.py reference) |
| ANA-04 | IASLC grading distribution across G1–G4 | IASLC grade NOT in adata.obs or ggo_config.yml — must derive from PATH STAGE column in clinical Excel or MM Class; requires patient-level merge with clinical metadata |
| ANA-05 | Count of PNS/IAC cases from Fig 1b that fall in G4 | Direct cross-tabulation from adata.obs: radio='PNS', pathology='IAC' → 8 patients; in adata 13 PNS in G4 (all PNS), need to filter to PNS-IAC subgroup |
| ANA-06 | Adjusted p-values for Fig 2l–m M1 macrophage N vs AIS | Requires loading the macrophage polarization density object from myeloid_analysis.py; imc.tl.grouped_mwu_test result; confirm BH correction was applied |
| ANA-07 | Inter-panel consistency between PANEL_G and PANEL_H documented | Correlate celltype densities between panels for shared cell types; document in analysis record (text note, not new figure) |
</phase_requirements>

---

## Summary

Phase 1 implements seven targeted computational analyses against the 122-patient GGO-IMC cohort, all starting from existing AnnData objects already on disk. The primary data object (`results/patient_celltype_broad_density_group.h5ad`, 122 patients × 14 cell-type densities) already contains the Group labels, mutation columns, radio, and pathology metadata needed for CLUS-01 through CLUS-03 and ANA-02 through ANA-05.

The largest technical risk is ANA-01 (Yoffe stratification): the Yoffe cohort AnnData is not referenced in `ggo_config.yml` and its local path is unknown. This must be discovered as Wave 0 work before centroid-based stratification can proceed. The process-score computation itself is well-understood — it mirrors the existing `sc.tl.score_genes` calls in `patient_group.py`, and the Yoffe dataset reportedly uses the same cell-type gene lists.

ANA-04 (IASLC grading) is a second moderate risk: the clinical Excel contains `PATH STAGE` (inconsistently formatted, 24+ variants) and `MM Class` (ais/mia/invasive) — neither maps directly to IASLC 2020 grade 1/2/3. A clean ordinal mapping rule must be established before plotting. The remaining analyses (ANA-05, ANA-06, ANA-07) are low-risk desk exercises: direct cross-tabulations, retrieving prior test results, and a short documentation note.

**Primary recommendation:** Write one new script `scripts/clustering_validation.py` for CLUS-01/02/03, one script `scripts/yoffe_stratification.py` for ANA-01, and one script `scripts/clinical_covariate_breakdown.py` for ANA-02–05. ANA-06 is a diagnostic check inside myeloid_analysis.py. ANA-07 is a prose note in an analysis record file.

---

## Standard Stack

### Core
| Library | Version | Purpose | Why Standard |
|---------|---------|---------|--------------|
| scanpy | 1.11.5 | Load/manipulate AnnData objects | Project-wide standard; already used in all scripts |
| numpy | 2.0.2 | Numerical operations, array math | Project-wide standard |
| pandas | (project) | Clinical metadata merging, cross-tabs | Project-wide standard |
| scipy.stats | 1.15.3 | chi2_contingency, fisher_exact, wilson CI | Already used in project; standard for these tests |
| sklearn.cluster | (scipy stack) | KMeans clustering | Standard; no hand-rolling |
| sklearn.metrics | (scipy stack) | adjusted_rand_score, silhouette_score | Standard; no hand-rolling |
| statsmodels.stats.proportion | (project scipy) | proportion_confint for Wilson/Clopper-Pearson CIs | Standard for CI computation; alternative to hand-rolling |
| matplotlib | 3.10.7 | Figure output | Project-wide standard |
| seaborn | 0.13.2 | Statistical bar plots | Project-wide standard |

### Supporting
| Library | Version | Purpose | When to Use |
|---------|---------|---------|-------------|
| scipy.spatial.distance | (scipy) | cdist for centroid-nearest-neighbor assignment | ANA-01 Yoffe stratification |
| imc_analysis | (custom) | grouped_mwu_test, celltype_density | ANA-06: retrieve existing test results |
| anndata | 0.11.4 | concat for multi-panel operations | ANA-07 inter-panel consistency |

**Installation:** All packages are already installed in the project environment. No new installations needed.

---

## Architecture Patterns

### Recommended Script Structure for New Analyses

New scripts follow the established pattern from `patient_group.py`:

```
scripts/
├── clustering_validation.py    # CLUS-01, CLUS-02, CLUS-03
├── yoffe_stratification.py     # ANA-01
├── clinical_covariate_breakdown.py  # ANA-02, ANA-03, ANA-04, ANA-05
└── [inline check in myeloid_analysis.py or separate diagnostic]  # ANA-06
```

Analysis records (non-figure outputs) go to:
```
figures/supplementary/
├── silhouette_analysis.pdf     # CLUS-01 output
└── covariate_<name>.pdf        # ANA-02–04 outputs
```

And a prose analysis record:
```
analysis_records/
├── clustering_validation.md    # ANA-07 inter-panel note, CLUS-01/02/03 numbers
```

### Pattern 1: Standard Script Skeleton
**What:** All new scripts load config, load data, compute, and save via utils.
**When to use:** Every new script in this phase.

```python
"""Module docstring describing purpose."""
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns
import matplotlib.pyplot as plt

from utils import load_config, save_figure, ensure_dir


if __name__ == '__main__':
    metadata = load_config()
    ensure_dir('figures/supplementary')

    adata = sc.read(
        metadata['patient_celltype_broad_clustered'],
        backup_url=metadata['patient_group_url'],
    )
    # ... analysis ...
    save_figure('figures/supplementary/silhouette_analysis.pdf', fig=fig)
```

### Pattern 2: Data Loading for Clustering AnnData
**What:** The primary patient-level AnnData (`patient_celltype_broad_density_group.h5ad`) is the ground truth.
**When to use:** CLUS-01, CLUS-02, CLUS-03, ANA-02, ANA-03, ANA-05.

```python
# Source: verified from patient_group.py (scripts/patient_group.py)
adata = sc.read(
    metadata['patient_celltype_broad_clustered'],
    backup_url=metadata['patient_group_url'],
)
# adata.obs already contains: Group, radio, pathology, EGFR, KRAS, RBM10, TP53, Smoking Status
# adata.X = cell-type densities (122 patients × 14 cell types)
# adata.var.index = ['B', 'CD8 T', 'CD4 T', 'T reg', 'Epi.-like', 'Epi. Prol.',
#                    'Mesen.-like', 'Fib.', 'Endo.', 'Mast', 'Mono.', 'Mac.', 'NK', 'PMN-MDSC']
```

### Pattern 3: Silhouette Analysis Over k Range
**What:** Compute mean silhouette width for each k; plot line; mark k=4 peak.
**When to use:** CLUS-01.

```python
# Source: sklearn docs; euclidean recommended for this dimensionality
from sklearn.metrics import silhouette_score
from sklearn.cluster import AgglomerativeClustering, KMeans
import numpy as np

sc.pp.scale(adata)
sc.pp.pca(adata)
X_pca = adata.obsm['X_pca'][:, :10]  # Use top-10 PCs, consistent with patient_group.py

k_range = range(2, 7)  # k=2..6
sil_scores = []
for k in k_range:
    labels = AgglomerativeClustering(n_clusters=k, linkage='ward').fit_predict(X_pca)
    sil_scores.append(silhouette_score(X_pca, labels))
```

**Distance metric decision (Claude's discretion):** Use euclidean distance on PCA-transformed data (consistent with how hierarchical clustering was originally run). The 14-feature space (cell-type densities) is not high-dimensional enough to warrant cosine distance.

### Pattern 4: K-means Concordance (CLUS-02)
**What:** Fit KMeans(k=4, n_init=10, random_state=0); compare to hierarchical Group labels via ARI and % agreement.
**When to use:** CLUS-02.

```python
from sklearn.cluster import KMeans
from sklearn.metrics import adjusted_rand_score

kmeans = KMeans(n_clusters=4, n_init=10, random_state=0)
kmeans_labels = kmeans.fit_predict(X_pca)

# Strip 'Group' prefix for numeric comparison
hier_labels = adata.obs['Group'].str.replace('Group', '').astype(int) - 1

ari = adjusted_rand_score(hier_labels, kmeans_labels)
# % agreement requires optimal label alignment (Hungarian algorithm or itertools.permutations for k=4)
from scipy.optimize import linear_sum_assignment
from sklearn.metrics import confusion_matrix
cm = confusion_matrix(hier_labels, kmeans_labels)
row_ind, col_ind = linear_sum_assignment(-cm)
pct_agreement = cm[row_ind, col_ind].sum() / len(hier_labels)
```

### Pattern 5: G4 Sensitivity/Specificity (CLUS-03)
**What:** Treat G4 membership as a binary predictor of PNS/IAC histology. Compute sens/spec with Wilson 95% CI.
**When to use:** CLUS-03.

**Critical clarification from data inspection:** The CONTEXT.md note about "n=8 PNS cases" refers specifically to PNS patients who are also histologically IAC (the key misclassification subgroup). Total PNS = 61 patients (61/122). The sens/spec computation uses this IAC-PNS group as the positive class.

```python
# Confirmed from data: adata.obs['radio']=='PNS' & adata.obs['pathology']=='IAC' => 8 patients
# Confirmed: adata.obs['Group']=='Group4' counts => 25 patients total in G4
positive_class = (adata.obs['radio'] == 'PNS') & (adata.obs['pathology'] == 'IAC')
in_g4 = adata.obs['Group'].str.replace('Group', '') == '4'

tp = (positive_class & in_g4).sum()
fn = (positive_class & ~in_g4).sum()
tn = (~positive_class & ~in_g4).sum()
fp = (~positive_class & in_g4).sum()

sensitivity = tp / (tp + fn)   # PNS-IAC captured in G4 / all PNS-IAC
specificity = tn / (tn + fp)   # non-PNS-IAC not in G4 / all non-PNS-IAC

from statsmodels.stats.proportion import proportion_confint
sens_ci = proportion_confint(tp, tp + fn, method='wilson')
spec_ci = proportion_confint(tn, tn + fp, method='wilson')
```

### Pattern 6: Clinical Covariate Bar Plots (ANA-02, ANA-03, ANA-04)
**What:** Stacked bar plots showing % per group, with chi-squared or Fisher's exact p-value.
**When to use:** ANA-02 (mutations), ANA-03 (smoking), ANA-04 (IASLC).

```python
from scipy.stats import chi2_contingency, fisher_exact

# Cross-tab Group vs covariate
ct = pd.crosstab(adata.obs['Group'], adata.obs['EGFR'])

# Fisher's exact if any cell < 5, else chi-squared
if (ct < 5).any().any():
    _, pval = fisher_exact(ct.values)
else:
    _, pval, _, _ = chi2_contingency(ct.values)

# Stacked bar (% per group)
ct_pct = ct.div(ct.sum(axis=1), axis=0) * 100
ct_pct.plot(kind='bar', stacked=True, ax=ax)
ax.set_title(f'p = {pval:.3g}')
```

### Pattern 7: Centroid-Based Yoffe Stratification (ANA-01)
**What:** Compute process-score centroids from IMC patients, assign each Yoffe patient to nearest centroid.
**When to use:** ANA-01.

```python
from scipy.spatial.distance import cdist

# Step 1: Compute process scores for IMC patients (same logic as patient_group.py)
processes = {
    'Epithelial score': ['Epi.-like', 'Epi. Prol.'],   # adjust per panel
    'Mesenchymal score': ['Mesen.-like'],
    'Fibrosis score': ['Fib.'],
    'Angiogenesis score': ['Endo.'],
    'B&T score': ['B', 'CD4 T', 'CD8 T'],
    'Pan-Immune score': ['CD4 T', 'CD8 T', 'Mac.', 'NK', 'T reg', 'Mast', 'Mono.', 'PMN-MDSC'],
}
# EMT score = Mesenchymal - Epithelial
for score_name, cell_types in processes.items():
    sc.tl.score_genes(adata, gene_list=cell_types, score_name=score_name, ctrl_size=20)
adata.obs['EMT score'] = adata.obs['Mesenchymal score'] - adata.obs['Epithelial score']

score_cols = ['EMT score', 'Epithelial score', 'Mesenchymal score',
              'Fibrosis score', 'Angiogenesis score', 'B&T score', 'Pan-Immune score']

# Step 2: Compute centroids per Group
group_labels = adata.obs['Group'].str.replace('Group', '').astype(int)
score_matrix = adata.obs[score_cols].values
centroids = {}
for g in [1, 2, 3, 4]:
    mask = (group_labels == g)
    centroids[g] = score_matrix[mask].mean(axis=0)
centroid_matrix = np.array([centroids[g] for g in [1, 2, 3, 4]])

# Step 3: Assign Yoffe patients (after computing same scores on Yoffe AnnData)
yoffe_scores = yoffe_adata.obs[score_cols].values  # same score columns
dists = cdist(yoffe_scores, centroid_matrix, metric='euclidean')
yoffe_groups = np.argmin(dists, axis=1) + 1  # 1-indexed groups
yoffe_adata.obs['Assigned_Group'] = yoffe_groups
```

### Anti-Patterns to Avoid

- **Don't re-implement silhouette from scratch:** `sklearn.metrics.silhouette_score` handles the per-sample computation; no hand-rolling needed.
- **Don't use Leiden clustering:** Locked out. Use only k-means for alternative clustering (CLUS-02).
- **Don't assume Group labels are already stripped of 'Group' prefix:** In adata.obs they appear as 'Group1', 'Group2', 'Group3', 'Group4' — strip before numeric comparison.
- **Don't assume process scores are pre-computed in adata.obs:** Verified — adata.obs has NO score columns. Scores must be computed via sc.tl.score_genes each run, or compute once and cache in new columns.
- **Don't assume IASLC grading is in adata.obs:** It is not. Must derive from PATH STAGE or MM Class column of the clinical Excel. PATH STAGE has inconsistent formatting (24+ variants). Plan a normalization step.
- **Don't re-implement Fisher's exact for 2×k tables:** Use `scipy.stats.fisher_exact` for 2×2, or `scipy.stats.chi2_contingency` for larger tables; never custom implementations.
- **Don't hard-code figure paths:** Use `ensure_dir()` and `save_figure()` from utils.py.

---

## Don't Hand-Roll

| Problem | Don't Build | Use Instead | Why |
|---------|-------------|-------------|-----|
| Silhouette width | Custom distance/cluster width loop | sklearn.metrics.silhouette_score | Handles edge cases; vectorized |
| ARI calculation | Custom label permutation score | sklearn.metrics.adjusted_rand_score | Chance-corrected standard; 1 line |
| Label alignment for concordance | Brute-force 24 permutations | scipy.optimize.linear_sum_assignment | Optimal O(n³) Hungarian assignment |
| Wilson/Clopper-Pearson CI | Custom formula | statsmodels.stats.proportion.proportion_confint | Exact implementation; method='wilson' or 'beta' |
| Fisher's exact test | Custom hypergeometric | scipy.stats.fisher_exact | Standard; handles 2×2 |
| Chi-squared test | Custom formula | scipy.stats.chi2_contingency | Returns statistic, p, dof, expected |
| Centroid nearest-neighbor | Loop with manual dist calc | scipy.spatial.distance.cdist | Vectorized; all pairs at once |

---

## Common Pitfalls

### Pitfall 1: Group Label Format Mismatch
**What goes wrong:** adata.obs['Group'] contains strings 'Group1'–'Group4'; k-means produces integer labels 0–3. Direct comparison fails silently or produces wrong ARI.
**Why it happens:** Two different labeling conventions exist in the same analysis.
**How to avoid:** Always strip 'Group' prefix and adjust to 0-indexed integers before passing to ARI or label comparison. Use `adata.obs['Group'].str.replace('Group', '').astype(int) - 1`.
**Warning signs:** ARI = 0 or very low even though clusters look visually similar.

### Pitfall 2: PNS Count Confusion
**What goes wrong:** ANA-05 / CLUS-03 refers to "8 PNS cases" which means PNS-radio AND IAC-pathology patients. But there are 61 total PNS patients. Using 61 as denominator gives wrong sensitivity.
**Why it happens:** "PNS cases" in the reviewer response means the radiologically misclassified subset (PNS radiology but invasive IAC histology — the diagnostic blind spot).
**How to avoid:** Filter `(radio=='PNS') & (pathology=='IAC')` for the true n=8 positive class. Verified from live data: 8 patients.
**Warning signs:** Suspiciously high sensitivity when n of positives is large.

### Pitfall 3: IASLC Grade Not in AnnData
**What goes wrong:** Script crashes with KeyError or produces NaN plots for ANA-04 because 'IASLC grade' or 'grade' is not in adata.obs.
**Why it happens:** IASLC grading is not available as a pre-computed column in the patient AnnData.
**How to avoid:** Load clinical Excel in ANA-04 script; derive ordinal grade from PATH STAGE or MM Class after normalization. Then merge by patient ID (GGO ID is the common key in adata.obs).
**Warning signs:** Missing column in adata.obs; adata.obs has only: pathology, radio, EGFR, KRAS, RBM10, TP53, Smoking Status, Gender, Race, ROI_area, Group.

### Pitfall 4: Yoffe AnnData Path Unknown
**What goes wrong:** ANA-01 script crashes because Yoffe data path is not in ggo_config.yml.
**Why it happens:** PROJECT.md says Yoffe is "already integrated" but the YAML config has no yoffe-related key (verified: grepped ggo_config.yml, no matches for 'yoffe', 'external', or 'cohort').
**How to avoid:** Wave 0 task must locate the Yoffe AnnData on disk (likely in `results/` or a separate directory) and add its path to ggo_config.yml before the stratification script is written.
**Warning signs:** ggo_config.yml has no Yoffe key.

### Pitfall 5: Process Score Computation Requires Matching Feature Names
**What goes wrong:** sc.tl.score_genes uses var.index (gene names / marker names) to look up the cell types in processes dict. If the Yoffe AnnData uses different feature names, the score computation fails or produces NaN.
**Why it happens:** The IMC AnnData uses cell-type names as var.index (e.g., 'Epi.-like', 'Fib.'), while Yoffe spatial transcriptomics might use gene names. The "same cell-type gene lists" in CONTEXT.md means Yoffe must be pre-processed to cell-type abundance space first.
**How to avoid:** Confirm Yoffe AnnData's var.index format before designing the score computation. If Yoffe uses gene-level data, use the actual marker gene lists from processes dict; if already cell-type abundance, use cell-type names directly.
**Warning signs:** sc.tl.score_genes raises "gene not found" warnings or returns uniform scores.

### Pitfall 6: PATH STAGE Inconsistency for IASLC Grade
**What goes wrong:** ANA-04 produces many NaN groups because PATH STAGE has 24+ format variants ('IA', 'taisn0', 't1an0', '1A', 'Stage IA', etc.) that aren't mapped to grade.
**Why it happens:** Clinical Excel data entry was inconsistent.
**How to avoid:** Use MM Class column as a simpler proxy (ais/mia/invasive — already clean with only 3 variants); or build explicit normalization dict for PATH STAGE. Prefer MM Class for plotting since it closely mirrors pathology categories already used. Note IASLC 2020 grading (1/2/3) technically requires additional histological data not available in this dataset — use PATH STAGE or MM Class as the closest surrogate.
**Warning signs:** Large number of NaN/unknown values when grouping by grade.

### Pitfall 7: Adjacent Multiple Testing for ANA-06
**What goes wrong:** Reporting raw p-values instead of BH-adjusted p-values for Fig 2l–m M1 comparisons, or confirming the wrong comparison.
**Why it happens:** myeloid_analysis.py calls imc.tl.grouped_mwu_test which may or may not apply multiple testing correction internally.
**How to avoid:** Check what imc.tl.grouped_mwu_test returns — specifically whether it stores raw or adjusted p-values. Apply `statsmodels.stats.multitest.multipletests` with method='fdr_bh' if correction was not already applied. The comparison of interest is M1 macrophage density: Normal tissue vs AIS.
**Warning signs:** p-values suspiciously uniform or suspiciously all below threshold.

---

## Code Examples

### Silhouette Loop with Line Plot
```python
# Source: sklearn.metrics.silhouette_score + patient_group.py data loading pattern
import matplotlib.pyplot as plt
from sklearn.metrics import silhouette_score
from sklearn.cluster import AgglomerativeClustering

sc.pp.scale(adata)
sc.pp.pca(adata)
X_pca = adata.obsm['X_pca'][:, :10]

k_range = range(2, 7)
sil_scores = []
for k in k_range:
    labels = AgglomerativeClustering(n_clusters=k, linkage='ward').fit_predict(X_pca)
    sil_scores.append(silhouette_score(X_pca, labels, metric='euclidean'))

fig, ax = plt.subplots(figsize=(4, 3), dpi=300)
ax.plot(list(k_range), sil_scores, marker='o', color='#2A363B')
ax.axvline(x=4, color='#E84A5F', linestyle='--', alpha=0.7, label='k=4 (selected)')
ax.set_xlabel('Number of clusters (k)')
ax.set_ylabel('Mean silhouette width')
ax.set_xticks(list(k_range))
ax.legend()
sns.despine()
```

### Stacked Bar with Fisher/Chi2 Annotation
```python
# Source: clinical_features.py reference + scipy.stats
from scipy.stats import chi2_contingency, fisher_exact
import matplotlib.pyplot as plt
import seaborn as sns

def covariate_bar_plot(adata_obs, group_col, covariate_col, ax, title):
    ct = pd.crosstab(adata_obs[group_col], adata_obs[covariate_col])
    if (ct < 5).values.any():
        # Use chi2 for multi-category (Fisher's only valid for 2x2)
        _, pval, _, _ = chi2_contingency(ct.values)
    else:
        _, pval, _, _ = chi2_contingency(ct.values)
    ct_pct = ct.div(ct.sum(axis=1), axis=0) * 100
    ct_pct.plot(kind='bar', stacked=True, ax=ax, legend=True)
    ax.set_title(f'{title}\np = {pval:.3g}', fontsize=9)
    ax.set_xlabel('Group')
    ax.set_ylabel('Proportion (%)')
    ax.tick_params(axis='x', rotation=0)
```

Note: Fisher's exact is strictly for 2×2 tables. For mutation data (WT/MUT × 4 groups), use chi2_contingency.

### Wilson CI for Proportion
```python
# Source: statsmodels docs
from statsmodels.stats.proportion import proportion_confint

n_success = tp   # e.g., PNS-IAC patients in G4
n_trials = tp + fn  # total PNS-IAC patients

ci_low, ci_high = proportion_confint(n_success, n_trials, method='wilson', alpha=0.05)
print(f"Sensitivity: {n_success/n_trials:.1%} (95% CI: {ci_low:.1%}–{ci_high:.1%})")
```

---

## State of the Art

| Old Approach | Current Approach | Impact |
|--------------|------------------|--------|
| Elbow method for k selection | Silhouette width as primary criterion | Silhouette is quantitative and directly reviewable; elbow is subjective |
| Single-method clustering | ARI concordance across methods | ARI is the community standard for cluster agreement in scRNA-seq papers |
| Point estimates only | CI-accompanied proportions | Required by reviewers; Wilson CI is well-behaved at small n |

---

## Open Questions

1. **Yoffe AnnData path and feature space**
   - What we know: PROJECT.md says it is "already integrated"; ggo_config.yml has no Yoffe key; Yoffe et al. is a spatial transcriptomics cohort (not IMC)
   - What's unclear: Whether Yoffe is stored as cell-type abundances (compatible with process score approach) or as raw gene expression (requires cell-type deconvolution step first)
   - Recommendation: Wave 0 task — locate Yoffe h5ad on disk (search `results/` and any external data dirs), inspect adata.var to confirm feature space, add path to ggo_config.yml

2. **IASLC grading definition for ANA-04**
   - What we know: Clinical Excel has PATH STAGE (24+ format variants) and MM Class (ais/mia/invasive, 49 patients with data out of 146 rows); IASLC 2020 grading (1/2/3) requires mitotic count + tumor necrosis not available here
   - What's unclear: Whether reviewers expect true IASLC 2020 grade or accept an ordinal proxy
   - Recommendation: Use MM Class (ais=non-invasive → grade surrogate, mia=minimally invasive, invasive=IAC) as the cleanest available grading proxy; flag in analysis notes that true IASLC grade is unavailable

3. **imc.tl.grouped_mwu_test multiple testing correction**
   - What we know: myeloid_analysis.py calls this function; result feeds Fig 2l-m; imc_analysis is a custom internal package
   - What's unclear: Whether BH correction is applied internally or must be applied post-hoc
   - Recommendation: ANA-06 task must inspect the test result object and check for a raw p-value column; if present, apply BH correction explicitly

4. **Number of PCs to use for clustering validation**
   - What we know: patient_group.py calls sc.pp.pca without specifying n_comps; default is 50 PCs for scanpy; X has only 14 features so PCA will produce max 13 meaningful components
   - What's unclear: How many PCs the original hierarchical clustering used
   - Recommendation: Use top-10 PCs for silhouette/k-means to match typical usage; document this choice in the analysis record

---

## Validation Architecture

### Test Framework
| Property | Value |
|----------|-------|
| Framework | pytest (project uses pytest; see TESTING.md) |
| Config file | `pyproject.toml` (`[tool.pytest.ini_options]`, `pythonpath = ["scripts"]`) |
| Quick run command | `pytest tests/test_imports.py -v` |
| Full suite command | `pytest -v` |

### Phase Requirements → Test Map
| Req ID | Behavior | Test Type | Automated Command | File Exists? |
|--------|----------|-----------|-------------------|-------------|
| CLUS-01 | Silhouette scores are numeric, len=5, max near k=4 | unit | `pytest tests/test_units.py::TestSilhouetteAnalysis -x` | ❌ Wave 0 |
| CLUS-02 | ARI in [0,1]; pct_agreement in [0,1] | unit | `pytest tests/test_units.py::TestKmeansConcordance -x` | ❌ Wave 0 |
| CLUS-03 | Sensitivity + specificity each in [0,1]; CI width > 0 | unit | `pytest tests/test_units.py::TestG4SensSpec -x` | ❌ Wave 0 |
| ANA-01 | All Yoffe patients assigned to group 1–4; chi2/Fisher p computable | smoke | `pytest tests/test_imports.py::test_import_yoffe_stratification -x` | ❌ Wave 0 |
| ANA-02 | Cross-tab shape correct; p-value is float | unit | `pytest tests/test_units.py::TestCovariateBreakdown -x` | ❌ Wave 0 |
| ANA-03 | Smoking covariate bar plot saved to disk | smoke | `pytest tests/test_imports.py::test_import_clinical_covariate_breakdown -x` | ❌ Wave 0 |
| ANA-04 | IASLC/grade column mapped from clinical Excel without NaN > 50% | unit | `pytest tests/test_units.py::TestIASLCMapping -x` | ❌ Wave 0 |
| ANA-05 | PNS-IAC count = 8 from live data (data assertion) | unit | `pytest tests/test_units.py::TestPNSIACCount -x` | ❌ Wave 0 |
| ANA-06 | Adjusted p-values are <= raw p-values or equal | manual-only | N/A — diagnostic inspection | N/A |
| ANA-07 | Analysis record file exists and contains inter-panel note | smoke | `pytest tests/test_imports.py::test_analysis_record_exists -x` | ❌ Wave 0 |

Note: The project's test suite currently mocks all data loading, so unit tests for new scripts can follow the same conftest.py mock pattern. Tests validate computation logic on fake data; actual data assertions (like ANA-05 n=8 check) are manual-only or integrated into the script as an assertion with a clear error message.

### Sampling Rate
- **Per task commit:** `pytest tests/test_imports.py -v`
- **Per wave merge:** `pytest -v`
- **Phase gate:** Full suite green before `/gsd:verify-work`

### Wave 0 Gaps
- [ ] `tests/test_units.py` — add test classes: TestSilhouetteAnalysis, TestKmeansConcordance, TestG4SensSpec, TestCovariateBreakdown, TestIASLCMapping, TestPNSIACCount
- [ ] `tests/test_imports.py` — add smoke tests for: clustering_validation, yoffe_stratification, clinical_covariate_breakdown
- [ ] `scripts/clustering_validation.py` — new file; covers CLUS-01, CLUS-02, CLUS-03
- [ ] `scripts/yoffe_stratification.py` — new file; covers ANA-01 (after Yoffe data path confirmed)
- [ ] `scripts/clinical_covariate_breakdown.py` — new file; covers ANA-02, ANA-03, ANA-04, ANA-05
- [ ] `analysis_records/` directory — new; hosts ANA-07 prose note and clustering validation numbers
- [ ] Yoffe AnnData path in ggo_config.yml — must be added before ANA-01 implementation

---

## Data Reality: Verified Facts from Live Inspection

These facts were verified by loading `results/patient_celltype_broad_density_group.h5ad` directly:

| Fact | Value |
|------|-------|
| Total patients | 122 |
| Group distribution | G1: 40, G2: 7, G3: 50, G4: 25 |
| Radio distribution | PNS: 61, PS: 31, S: 24 |
| Pathology distribution | AIS: 44, MIA: 43, IAC: 35 |
| PNS+IAC patients (key diagnostic blind spot) | 8 |
| PNS patients in G4 | 13 (out of 25 total G4; out of 61 total PNS) |
| Mutation columns in adata.obs | EGFR, KRAS, RBM10, TP53 (as 'WT'/'MUT' Categorical) |
| Process scores pre-computed | NO — must compute via sc.tl.score_genes each time |
| adata.X content | Cell-type densities (cells/ROI_area units) |
| adata.var.index | 14 broad cell types |

---

## Sources

### Primary (HIGH confidence)
- Direct code inspection: `scripts/patient_group.py`, `scripts/utils.py`, `scripts/myeloid_analysis.py`, `scripts/exploratory/mutation.py`, `scripts/exploratory/clinical_features.py`, `scripts/exploratory/create_patient_density_matrix.py`
- Direct data inspection: `results/patient_celltype_broad_density_group.h5ad` (loaded and queried)
- Direct metadata inspection: `metadata/ggo_config.yml`, `metadata/De-identified J&J Clinical Annotations_UpdatedAAH031523 (1).xlsx`
- `.planning/codebase/STACK.md`, `.planning/codebase/CONVENTIONS.md`, `.planning/codebase/ARCHITECTURE.md`, `.planning/codebase/TESTING.md`, `.planning/codebase/CONCERNS.md`
- `metadata/mutation_piv.csv` — confirmed EGFR, KRAS, RBM10, TP53 columns present

### Secondary (MEDIUM confidence)
- sklearn documentation for silhouette_score, adjusted_rand_score, KMeans — standard library behavior
- statsmodels proportion_confint — Wilson/Clopper-Pearson CI; well-established

### Tertiary (LOW confidence)
- Yoffe cohort feature space — inferred from CONTEXT.md ("same cell-type gene lists"); unverified until data located

---

## Metadata

**Confidence breakdown:**
- Standard stack: HIGH — all packages already in project; versions verified from STACK.md
- Architecture / script structure: HIGH — directly mirrors existing pattern in patient_group.py
- Data facts (counts, columns): HIGH — verified by loading live h5ad
- Pitfalls: HIGH — identified from direct code and data inspection
- Yoffe data path: LOW — not found in config; must be discovered

**Research date:** 2026-03-20
**Valid until:** 2026-06-20 (stable domain; data is static)
