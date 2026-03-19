# ggo-imc Manuscript Revision

## What This Is

Revision of a published Cancer Cell (2025) paper: "Spatial profiling of early-stage lung adenocarcinoma reveals patterns of immunomodulation and epithelial plasticity." The paper performed IMC analysis on 122 early-stage LUAD patients, characterizing immune, stromal, and epithelial trajectories across AAH → AIS → MIA → IAC progression, and identifying 4 patient groups (G1–G4) with distinct clinical phenotypes. The revision addresses computational reviewer requests and, in a second milestone, pursues a broader design and storytelling overhaul to strengthen the manuscript's positioning.

## Core Value

Reviewers want rigor and novelty. The resubmission must demonstrate analytical robustness (patient clustering, statistics), clearly distinguish findings from Zhu et al. 2025, and fix all figure clarity issues — without new wet-lab experiments.

## Requirements

### Validated

- ✓ IMC data loaded and QC-filtered — existing
- ✓ Cell types assigned (35–37 marker panels, immune + tumor) — existing
- ✓ All 5 manuscript figures generated — existing
- ✓ Spatial domain (UTAG) analysis complete — existing
- ✓ EMP validation via Yoffe et al. spatial transcriptomics — existing

### Active

#### Milestone 1 — Reviewer Response (Resubmission)

**Clustering & Robustness (Fig 5)**
- [ ] Add silhouette width analysis to justify k=4 patient grouping
- [ ] Apply alternative clustering method (k-means or Leiden) and compare to hierarchical result
- [ ] Report sensitivity/specificity of G4 for identifying PNS/IAC misclassified cases

**New Analyses**
- [ ] Stratify Yoffe et al. cohort into G1–G4 using this study's marker signatures
- [ ] Driver mutation distribution across G1–G4 (EGFR, KRAS, RBM10, TP53)
- [ ] Smoking status distribution across G1–G4
- [ ] IASLC grading across groups, especially enrichment in G4
- [ ] PNS/IAC count breakdown among the 8 PNS/IAC cases shown in Fig 1b → how many fall in G4
- [ ] Clarify MDM representation within Fig 2L myeloid categories (text/methods fix)
- [ ] Confirm/report adjusted p-values for Fig 2l–m (M1 macrophage N vs AIS comparison)

**Figure Fixes**
- [ ] Fig 2g, 2h, 2m: add individual data points to bar plots
- [ ] Fig 4d/4e: clarify score scale meaning; restrict to lymphoid markers in TLS clusters, myeloid in macrophage clusters
- [ ] Fig 4f: fix label direction (tumour vs normal vs normal vs tumour)
- [ ] Fig 5d/e: add PNS/PS/S percentage breakdown per group (link to Table 2)
- [ ] Extended Fig 4: fix legend (two conflicting panel a/b legends; missing violin + clinical tumor size plots)
- [ ] Extended Fig 7: add missing panel c
- [ ] Standardize GranzymeB → GZMB nomenclature throughout

**Text & Writing**
- [ ] Replace all causal/mechanistic language with correlational language (all 3 reviewers flagged)
- [ ] Cite and discuss Zhu et al. 2025 (Cancer Cell, PMID 40345189) — overlap and differentiation
- [ ] Fix typo: 'Grou p 3' (line 263)
- [ ] Fix: line 294 references G3 but should reference G4 for misclassification
- [ ] Fix: line 231 'black and green dotted lines'
- [ ] Replace 'Clara cells' with 'Club cells' (line 37)
- [ ] Define PanCK and other key markers at first mention
- [ ] Expand figure legends throughout for clarity
- [ ] Add methods description for spatial spot program scoring (AT1/AT2/EMP)
- [ ] Fix line 135: "Similar to our previous findings" (referenced paper is not from authors)
- [ ] Address Fig 4f axis ambiguity in text and legend
- [ ] Remove/reframe US statistics anchor in introduction
- [ ] Clarify hallmark process basis in Fig 5c (without RNA-seq)
- [ ] Clarify inter-panel consistency between tumor and immune IMC panels (Rev2 #1)

#### Milestone 2 — Manuscript Overhaul (Future)

- [ ] Redesign figure aesthetics and layout for impact
- [ ] Restructure narrative arc to foreground novelty vs Zhu et al. 2025
- [ ] Sharpen abstract and title for positioning
- [ ] Strengthen Fig 5 as the key novel claim with cleaner visual story
- [ ] Consider panel reorganization to reduce figure complexity

### Out of Scope

- Additional wet-lab staining (FAP, PDGFRβ, CXCL12 for CAF subtypes) — requires new experiments
- Functional/mechanistic validation of EMP or immune exclusion — requires new experiments
- Held-out independent validation cohort — not available

## Context

- **Platform:** IMC (Imaging Mass Cytometry), TMA format
- **Cohort:** 122 early-stage LUAD patients; AAH, AIS, MIA, IAC + paired normal
- **Panels:** 35–37 antibodies each; tumor panel + immune panel on adjacent sections
- **Pipeline:** Snakemake + Python (sc_tools), AnnData objects
- **External validation dataset:** Yoffe et al. spatial transcriptomics (already integrated)
- **Reviewers:** 3 reviewers — R1 skeptical of novelty, R2 technically constructive, R3 most critical on robustness and causal overclaiming
- **Draft response:** `manuscript/Reviewer Response.docx` — partially written

## Constraints

- **No new experiments:** All improvements must be purely computational/analytical
- **Data:** Existing AnnData + metadata only; Yoffe et al. dataset already loaded
- **R dependency:** `scripts/asd.R` (Fig 5 patient analysis) runs outside container in `ggo_imc` conda env
- **Runtime:** `SC_TOOLS_RUNTIME=none snakemake ...` with `ggo_imc` conda env

## Key Decisions

| Decision | Rationale | Outcome |
|----------|-----------|---------|
| Two-milestone revision strategy | Decouple reviewer response (resubmission) from design overhaul (longer-term improvement) | — Pending |
| No new wet-lab experiments | Scope constraint from authors; focus on computational improvements only | — Pending |
| Milestone 1 priority: clustering robustness | R2 and R3 both flag Fig 5 patient grouping as weakest analytical point | — Pending |

---
*Last updated: 2026-03-19 after project initialization from reviewer comments*
