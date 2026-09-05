# Limitations

- Only **two Gate A PASS proteins** (CYP2C9 and PTEN) support the binary primary benchmark; GCK and VKORC1 are descriptive-only ([cohort counts](../results/metadata/task2_counts.csv); sample unit: variants/positions).
- The evidence spans four assay systems/readouts, with abundance and function measured in different experimental contexts ([assay inventory](../results/metadata/assay_inventory.csv); sample unit: assays).
- The primary endpoint conditions on measured abundance being preserved. It is not a sequence-only deployment model and does not estimate function in unmeasured proteins or variants.
- Label uncertainty differs by source: CYP2C9 and VKORC1 use author categories, PTEN uses a high-confidence score rule, and GCK uses score/SE confidence intervals ([analysis freeze](../config/analysis_freeze.yaml); sample unit: labeled variants).
- There is no clinical endpoint, and the analysis does not establish pathogenicity, causality, stability, a causal mechanism, or allostery.
- Confidence intervals are position bootstraps over measured residue positions; they do not cover unmeasured positions or all assay error ([bootstrap summary](../results/metadata/task4_bootstrap_summary.csv); sample unit: measured positions).
- LOPO has only two proteins, so held-out-protein estimates are a two-protein comparative check, not proteome generalization ([task 5 metrics](../results/metadata/task5_metrics.csv); sample unit: proteins/variants).
- Upstream cross-assay score discrepancies are real and retained as an audit/sensitivity source; function-assay scores were not silently coalesced with the abundance-assay production source ([score consistency](../results/metadata/score_consistency.csv); sample unit: protein/model groups).
- The v0.1 scope is the amended primary binary benchmark only; Spearman score/function and grouped-OOF isotonic-residual endpoints were deferred after the primary audit, so no secondary inference is reported ([cohort freeze](../config/cohort_freeze.json); sample units: variants, residue positions, and proteins).
