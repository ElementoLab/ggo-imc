---
phase: 1
slug: analytical-robustness
status: draft
nyquist_compliant: false
wave_0_complete: false
created: 2026-03-20
---

# Phase 1 — Validation Strategy

> Per-phase validation contract for feedback sampling during execution.

---

## Test Infrastructure

| Property | Value |
|----------|-------|
| **Framework** | pytest (present in sc_tools env) + file existence checks |
| **Config file** | none — Wave 0 stubs created per task |
| **Quick run command** | `python -c "import h5py; import anndata; print('ok')"` (smoke test) |
| **Full suite command** | `pytest tests/phase1/ -v` (Wave 0 creates) |
| **Estimated runtime** | ~30 seconds |

---

## Sampling Rate

- **After every task commit:** Run quick smoke test (file exists + output shape check)
- **After every plan wave:** Run full suite command
- **Before `/gsd:verify-work`:** Full suite must be green
- **Max feedback latency:** 60 seconds

---

## Per-Task Verification Map

| Task ID | Plan | Wave | Requirement | Test Type | Automated Command | File Exists | Status |
|---------|------|------|-------------|-----------|-------------------|-------------|--------|
| silhouette-plot | 01 | 1 | CLUS-01 | file check | `ls figures/figure5/silhouette_k_sweep.pdf` | ❌ W0 | ⬜ pending |
| kmeans-concordance | 01 | 1 | CLUS-02 | file check | `ls figures/figure5/kmeans_concordance.txt` | ❌ W0 | ⬜ pending |
| g4-sens-spec | 01 | 1 | CLUS-03 | file check | `ls results/g4_sensitivity_specificity.txt` | ❌ W0 | ⬜ pending |
| yoffe-stratification | 02 | 2 | ANA-01 | file check | `ls results/yoffe_group_assignments.csv` | ❌ W0 | ⬜ pending |
| mutation-breakdown | 03 | 2 | ANA-02 | file check | `ls figures/supplementary/mutation_by_group.pdf` | ❌ W0 | ⬜ pending |
| smoking-breakdown | 03 | 2 | ANA-03 | file check | `ls figures/supplementary/smoking_by_group.pdf` | ❌ W0 | ⬜ pending |
| iaslc-breakdown | 03 | 2 | ANA-04 | file check | `ls figures/supplementary/iaslc_by_group.pdf` | ❌ W0 | ⬜ pending |
| pns-iac-counts | 03 | 2 | ANA-05 | file check | `ls results/pns_iac_group_counts.txt` | ❌ W0 | ⬜ pending |
| adjusted-pvals | 04 | 2 | ANA-06 | file check | `ls results/m1_adjusted_pvals.txt` | ❌ W0 | ⬜ pending |
| interpanel-check | 04 | 2 | ANA-07 | manual | review analysis record | N/A | ⬜ pending |

*Status: ⬜ pending · ✅ green · ❌ red · ⚠️ flaky*

---

## Wave 0 Requirements

- [ ] Locate Yoffe AnnData on disk and add path to `metadata/ggo_config.yml` — required before ANA-01
- [ ] Confirm `imc.tl.grouped_mwu_test` BH correction behavior for ANA-06
- [ ] Verify `MM Class` column as IASLC proxy for ANA-04

*If none: "Existing infrastructure covers all phase requirements."*

---

## Manual-Only Verifications

| Behavior | Requirement | Why Manual | Test Instructions |
|----------|-------------|------------|-------------------|
| Inter-panel consistency analysis record | ANA-07 | Narrative/methods text, not a file output | Read analysis notes; confirm tumor vs immune panel consistency addressed in writing |
| IASLC proxy defensibility | ANA-04 | Judgment call on MM Class acceptability | Check if MM Class maps cleanly to IASLC 2020; document in analysis notes |

---

## Validation Sign-Off

- [ ] All tasks have `<automated>` verify or Wave 0 dependencies
- [ ] Sampling continuity: no 3 consecutive tasks without automated verify
- [ ] Wave 0 covers all MISSING references
- [ ] No watch-mode flags
- [ ] Feedback latency < 60s
- [ ] `nyquist_compliant: true` set in frontmatter

**Approval:** pending
