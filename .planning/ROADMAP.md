# Roadmap: ggo-imc Manuscript Revision

## Overview

This roadmap covers Milestone 1 of the ggo-imc revision: a computational response to three reviewers of a Cancer Cell paper on spatial IMC profiling of early-stage LUAD. The work groups naturally into three phases: new computations that produce results (Phase 1), figure panel revisions that incorporate those results (Phase 2), and manuscript text corrections that frame everything accurately (Phase 3). Phases must execute in this order — downstream figures depend on analytical outputs, and text depends on both.

## Phases

**Phase Numbering:**
- Integer phases (1, 2, 3): Planned milestone work
- Decimal phases (2.1, 2.2): Urgent insertions (marked with INSERTED)

Decimal phases appear between their surrounding integers in numeric order.

- [ ] **Phase 1: Analytical Robustness** - New computations: clustering validation, clinical subgroup breakdowns, and external cohort stratification
- [ ] **Phase 2: Figure Revisions** - All figure panel changes incorporating Phase 1 results plus visual fixes
- [ ] **Phase 3: Manuscript Text** - Writing corrections: causal language, citations, typos, legends, and methods

## Phase Details

### Phase 1: Analytical Robustness
**Goal**: All new computational analyses are complete, validated, and ready to feed into figures and text
**Depends on**: Nothing (first phase)
**Requirements**: CLUS-01, CLUS-02, CLUS-03, ANA-01, ANA-02, ANA-03, ANA-04, ANA-05, ANA-06, ANA-07
**Success Criteria** (what must be TRUE):
  1. Silhouette width plot exists confirming k=4 is optimal; alternative clustering (k-means or Leiden) concordance with hierarchical result is quantified and documented
  2. Sensitivity and specificity of G4 for identifying PNS/IAC cases are computed and a number can be cited in the response letter
  3. Yoffe et al. cohort is stratified into G1–G4 using this study's signatures; risk patterns in external cohort match internal findings directionally
  4. Clinical covariate tables (driver mutations, smoking, IASLC grade, PNS/IAC case counts) broken down by G1–G4 are produced and ready to insert into figures or supplementary
  5. Adjusted p-values for Fig 2l–m M1 comparisons are confirmed and inter-panel consistency between tumor and immune panels is addressed in the analysis record
**Plans**: TBD

### Phase 2: Figure Revisions
**Goal**: All figure panels in the manuscript correctly display individual data points, accurate labels, marker-appropriate scores, and complete extended figures
**Depends on**: Phase 1
**Requirements**: FIG-01, FIG-02, FIG-03, FIG-04, FIG-05, FIG-06, FIG-07
**Success Criteria** (what must be TRUE):
  1. Fig 2g, 2h, 2m show individual data points overlaid on bar plots (reviewer can verify from PDF)
  2. Fig 4d/4e score scale is labeled with a clear definition; TLS cluster display is restricted to lymphoid markers, macrophage cluster to myeloid markers
  3. Fig 4f axis labels correctly and consistently indicate tumour vs normal direction
  4. Fig 5d/5e include PNS/PS/S percentage breakdown per group with reference to Table 2
  5. Extended Fig 4 legend is corrected (no duplicate panel labels; missing violin and clinical tumor size plots added or legend updated to match); Extended Fig 7 includes panel c; GranzymeB is rendered as GZMB throughout all figures
**Plans**: TBD

### Phase 3: Manuscript Text
**Goal**: The manuscript text is factually corrected, uses appropriate correlational language, cites Zhu et al. 2025, and has complete figure legends and methods
**Depends on**: Phase 2
**Requirements**: TXT-01, TXT-02, TXT-03, TXT-04, TXT-05, TXT-06, TXT-07, TXT-08, TXT-09, TXT-10, TXT-11, TXT-12, TXT-13
**Success Criteria** (what must be TRUE):
  1. No sentence in the manuscript asserts causation or mechanism; all such language is replaced with correlational phrasing ("associated with", "correlated with") throughout
  2. Zhu et al. 2025 (PMID 40345189) is cited in the introduction and/or discussion with explicit articulation of overlap and how this study is differentiated
  3. All discrete text errors are fixed: 'Grou p 3' (line 263), G3→G4 (line 294), 'black and green dotted lines' (line 231), 'Clara cells'→'Club cells' (line 37), "Similar to our previous findings" (line 135)
  4. Figure legends are expanded with sufficient detail for a reader to interpret each panel without consulting methods; PanCK and key markers are defined at first use; MDM representation in Fig 2L and hallmark basis in Fig 5c are explained
  5. Methods section includes spatial spot program scoring (AT1/AT2/EMP) description; introduction US statistics anchor is softened or removed; inter-panel consistency is addressed in methods/results prose
**Plans**: TBD

## Progress

**Execution Order:**
Phases execute in numeric order: 1 → 2 → 3

| Phase | Plans Complete | Status | Completed |
|-------|----------------|--------|-----------|
| 1. Analytical Robustness | 0/TBD | Not started | - |
| 2. Figure Revisions | 0/TBD | Not started | - |
| 3. Manuscript Text | 0/TBD | Not started | - |
