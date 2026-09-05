# Release checklist

Statuses are evidence-backed and apply to this repository snapshot.

| Check | Status | Evidence / reason |
|---|---|---|
| Scientific analysis complete for amended v0.1 scope | PASS | `config/cohort_freeze.json` records the dated amendment; primary metrics, learned metrics, sensitivity table, and six figures are tracked. Deferred secondary endpoints are not claimed. |
| Primary cohort counts and Gate A scope | PASS | `task2_counts.csv`: only CYP2C9/PTEN pass; GCK/VKORC1 are descriptive-only. |
| Direct score metrics and CIs | PASS | `task4_primary_metrics.csv` and `task4_bootstrap_summary.csv`; CIs are position bootstraps with invalid fractions recorded. |
| Structure/sensitivity outcomes | PASS | `task5_metrics.csv`, `task5_paired_bootstrap_summary.csv`, and `task5_sensitivity_results.csv`; NOT_APPLICABLE controls retained. |
| Cohort freeze and secondary-endpoint amendment | PASS | `config/cohort_freeze.json` records 2,055 primary variants, 550 measured residue positions, and the 2026-09-05 post-analysis deferral of Spearman and grouped-OOF isotonic-residual endpoints. |
| Secondary continuous endpoints | NOT_APPLICABLE | Explicitly deferred from amended v0.1 scope; no Spearman or grouped-OOF isotonic-residual inference is claimed. |
| Six figure manifest/readability | PASS | Six PNGs are tracked with exact dimensions and SHA256 in `results/metadata/figure_manifest.csv`; original-size inspection completed. |
| Offline reproducibility and CI definition | PASS | `tests/test_reproducibility.py` and `.github/workflows/ci.yml`; Windows/Ubuntu, Python 3.11, frozen uv lock. |
| Row-level source-derived outputs externalized | PASS | Local `D:\protein_assay_work\release_artifacts\manifest.json` records seven intact copies, hashes, bytes, and rows; exact public-tree ignores are committed. |
| Public redistribution of upstream raw/row data | NOT_APPLICABLE | Upstream data permissions are not established, so raw archives, PDBs, author tables, audit samples, case tables, and all variant-level outputs are excluded from the public history. |
| Public push | PASS | The portfolio release uses a reviewed clean history containing code, aggregate results, figures, and provenance only. |
| Clinical/causal interpretation | NOT_APPLICABLE | No clinical endpoint or causal-mechanism analysis is in scope. |
