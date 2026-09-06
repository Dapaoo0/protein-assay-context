# Protein Assay Context

> **When measured abundance looks normal, can a protein still lose function—and can sequence models detect it?**

## Abstract

**Background.** Variant-effect predictors are often evaluated against a single assay, even though protein abundance and protein function measure different biological consequences. A missense variant can preserve measured abundance yet impair a separate functional readout, creating a stringent test of whether a model captures functional constraint beyond gross abundance loss.

**Methods.** We paired abundance and function deep-mutational-scanning measurements for 18,207 variants across CYP2C9, PTEN, GCK, and VKORC1. A frozen eligibility gate retained proteins with adequate positive/negative variant and residue-position support. ProteinGym ESM1v and ESM2 zero-shot scores were evaluated on abundance-preserved variants using position-clustered bootstrap intervals. Learned comparisons used nested residue-position-grouped cross-validation; experimental-structure descriptors were added only inside leakage-controlled folds.

**Results.** CYP2C9 and PTEN passed the primary gate, yielding **2,055 variants across 550 measured residue positions**. ESM2 achieved average precision (AP) **0.838 versus prevalence 0.594** for CYP2C9 and **0.559 versus 0.156** for PTEN. ESM1v reached AP 0.794 and 0.543. Adding structure descriptors to the ESM1v score-only baseline changed within-protein AP by **−0.087 [−0.163, −0.012]** for CYP2C9 and **−0.021 [−0.116, 0.071]** for PTEN; held-out-protein intervals included zero.

**Interpretation.** Pretrained sequence scores contain useful signal for functional impairment that is not reducible to measured abundance loss. In this benchmark, compact structure descriptors did not provide consistent incremental value. The result is assay-specific evidence—not a clinical pathogenicity prediction, causal mechanism, stability claim, or proteome-wide conclusion.

## Study design

![Cohort construction from paired abundance and function assays](results/figures/figure1_cohort_flow.png)

*Figure 1. Four paired-assay cohorts contributed 18,207 variants. Frozen class and residue-position gates restricted confirmatory benchmarking to CYP2C9 and PTEN; GCK and VKORC1 remain descriptive.*

The primary endpoint is **function-impaired versus function-preserved among variants with preserved measured abundance**. Grouping every split and uncertainty calculation by residue position prevents substitutions at the same site from leaking across evaluation folds.

| Primary cohort | CYP2C9 | PTEN |
|---|---:|---:|
| Variants | 694 | 1,361 |
| Measured positions | 262 | 288 |
| Impaired / preserved | 412 / 282 | 212 / 1,149 |
| Prevalence | 0.594 | 0.156 |
| ESM1v AP | 0.794 | 0.543 |
| ESM2 AP | **0.838** | **0.559** |

## Main result: sequence models recover hidden functional loss

![Direct zero-shot benchmark for ESM1v and ESM2](results/figures/figure4_predictor_forest.png)

*Figure 4. Direct zero-shot AP with 95% position-clustered bootstrap intervals. Prevalence is the no-skill reference for each protein; intervals describe variation over measured residue positions, not independent variants or assay uncertainty.*

The improvement above prevalence is especially pronounced for PTEN because function-impaired variants are relatively rare in its abundance-preserved cohort. ESM2 exceeds ESM1v by **+0.0426 AP [0.0192, 0.0687]** for CYP2C9, while the PTEN contrast of **+0.0165 [−0.0408, 0.0730]** remains uncertain. The evidence therefore supports protein-specific performance rather than a universal ranking claim.

## Does structure add predictive value?

![Within-protein and held-out-protein generalization](results/figures/figure5_generalization.png)

*Figure 5. Leakage-controlled comparisons of sequence-only and structure-added feature sets. B2 and B4 use ESM1v; the direct ESM2 result above is a separate zero-shot benchmark.*

Experimental structures mapped to 92.45% of CYP2C9 and 76.18% of PTEN reference positions. Nevertheless, adding accessibility, secondary structure, contact, confidence, and functional-site-distance descriptors did not consistently improve AP. This negative result is scientifically useful: structural availability alone does not guarantee incremental predictive information, especially in a two-protein benchmark with incomplete residue coverage.

## Research safeguards

- Cohorts, labels, score orientation, model comparison, and eligibility gates are frozen and machine-checked.
- Model coverage is 100% for both primary proteins; missing structural residues remain missing.
- Nested grouped-position evaluation prevents residue-level leakage.
- Uncertainty uses 2,000-draw paired position bootstraps, preserving the measured-position sampling unit.
- All raw and variant-level source-derived tables remain outside the public Git history; tracked results are aggregates, figures, and provenance.

## Reproduce and audit

```powershell
uv sync --frozen
uv run python -m pytest -q
uv run ruff check src tests scripts
uv build
```

The offline test suite validates cohort contracts, canonical variant joins, grouped evaluation, deterministic bootstrap behavior, figure hashes, and public-package exclusions. Reproducing the full analysis requires separately downloaded, checksum-verified official inputs.

[Detailed results](docs/RESULTS.md) · [Reproducibility](docs/REPRODUCIBILITY.md) · [Data provenance](docs/DATA_PROVENANCE.md) · [Limitations](docs/LIMITATIONS.md) · [Release checklist](docs/RELEASE_CHECKLIST.md) · [Novelty audit](docs/NOVELTY_AUDIT.md)

## Scope

This repository is a reproducible, assay-conditioned comparative case study. Only CYP2C9 and PTEN support the primary binary benchmark; GCK and VKORC1 are descriptive-only. The analysis does not establish clinical impact, molecular causality, allostery, stability, or generalization beyond the assayed proteins and experimental systems.

## License

Code is released under the [MIT License](LICENSE). Upstream datasets and structures retain their original terms and are not redistributed here.
