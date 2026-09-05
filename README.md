# Protein Assay Context

## Result first

Pretrained sequence scores recover functional loss among variants with preserved measured abundance in the two proteins that passed the primary cohort gate: **2,055 primary variants across 550 measured residue positions** (CYP2C9: 694 variants/262 positions; PTEN: 1,361 variants/288 positions). Direct **ESM2** average precision (AP) was **0.838 versus prevalence 0.594 for CYP2C9**, and **0.559 versus 0.156 for PTEN** ([primary metrics](results/metadata/task4_primary_metrics.csv), [Figure 4](results/figures/figure4_predictor_forest.png)). AP is the area under the precision–recall curve; prevalence is the positive fraction in the abundance-preserved cohort.

Here, abundance means measured protein amount in an assay, while function means a separate experimental readout such as probe binding or complementation. A function-impaired label therefore does not by itself identify a stability defect, mechanism, clinical effect, or causality.

Structure features did not improve the score-only **B2** baseline in within-protein evaluation: adding structure (**B4**) changed AP by **−0.087 [−0.163, −0.012]** for CYP2C9 and **−0.021 [−0.116, 0.071]** for PTEN ([paired contrasts](results/metadata/task5_paired_bootstrap_summary.csv), [Figure 5](results/figures/figure5_generalization.png)). B2/B4 use **ESM1v**; the direct ESM2 result above is a separate zero-shot model and is not the B2/B4 score.

## Figures

[Cohort flow](results/figures/figure1_cohort_flow.png) · [Assay context](results/figures/figure2_abundance_function_scatter.png) · [Phenotype fractions](results/figures/figure3_phenotype_fractions.png) · [Direct benchmark](results/figures/figure4_predictor_forest.png) · [Generalization](results/figures/figure5_generalization.png) · [Structure and cases](results/figures/figure6_structure_error_map.png)

## Quick start

```powershell
uv sync --frozen
uv run python -m pytest -q
uv run python -m ruff check src tests scripts
uv build
```

These commands are offline after dependency setup. The full pipeline requires separately downloaded, checksum-verified official inputs on `D:\protein_assay_work`; see [Reproducibility](docs/REPRODUCIBILITY.md).

## Read next

- [Results](docs/RESULTS.md) — direct zero-shot, learned OOF/LOPO, structure, and sensitivity evidence.
- [Limitations](docs/LIMITATIONS.md) — scope and interpretation boundaries.
- [Data provenance](docs/DATA_PROVENANCE.md) — source revisions, hashes, roles, and redistribution policy.
- [Release checklist](docs/RELEASE_CHECKLIST.md) — evidence-backed readiness status.
- [Novelty audit](docs/NOVELTY_AUDIT.md) — prior work and the leakage-controlled phenotype-conditioned extension framing.

## Scope

This is an assay-specific comparative case study, not a clinical pathogenicity model, causal-mechanism study, allostery analysis, stability predictor, or proteome-wide claim. Only CYP2C9 and PTEN pass the binary primary gate; GCK and VKORC1 remain descriptive-only. Row-level source-derived tables, audit samples, case tables, and raw archives stay outside the public tree; compact aggregates, figures, and provenance are tracked.
