# Results

## Cohort and endpoint

The paired missense cohort contains 4,421 CYP2C9, 4,839 PTEN, 8,255 GCK, and 692 VKORC1 variants ([cohort counts](../results/metadata/task2_counts.csv); sample unit: variants). The binary endpoint is function-impaired versus function-preserved among variants with preserved measured abundance. CYP2C9 and PTEN pass Gate A; GCK and VKORC1 are descriptive-only because their class/position minima are not met. The primary benchmark therefore contains 694 CYP2C9 variants across 262 positions and 1,361 PTEN variants across 288 positions ([primary metrics](../results/metadata/task4_primary_metrics.csv); sample units: variants and residue positions).

The v0.1 release scope is frozen in [`config/cohort_freeze.json`](../config/cohort_freeze.json): the primary binary benchmark was completed first, and a dated post-analysis amendment defers Spearman score/function and grouped-OOF isotonic-residual endpoints. No secondary inference from those endpoints is claimed in this release.

## Direct zero-shot scores

The frozen direct comparison uses ProteinGym ESM1v and ESM2 zero-shot scores, oriented so larger values indicate impairment. ESM2 AP is 0.838 for CYP2C9 (prevalence 0.594; AP lift 1.411) and 0.559 for PTEN (prevalence 0.156; AP lift 3.587). ESM1v AP is 0.794 and 0.543, respectively ([task 4 primary metrics](../results/metadata/task4_primary_metrics.csv); sample unit: primary variants). The 95% intervals are **2,000-draw position bootstraps**: whole measured-variant clusters at a residue position are sampled together, so they quantify variation over measured positions rather than independent-variant or assay uncertainty. Per-protein heterogeneity is material: PTEN has lower AP but much lower prevalence, producing the larger AP lift; the equal-protein macro AP is 0.668 for ESM1v and 0.698 for ESM2 ([same table](../results/metadata/task4_primary_metrics.csv); sample unit: proteins for macro rows).

## Learned OOF and held-out-protein evaluation

Within-protein nested grouped-position OOF results show B2 (ESM1v score-only) AP of 0.785 for CYP2C9 and 0.498 for PTEN; B4 (ESM1v plus structure) is 0.697 and 0.477 ([learned metrics](../results/metadata/task5_metrics.csv); sample unit: variants, with position-grouped folds). The corresponding held-out-protein rows are CYP2C9 0.794/0.792 (B2/B4) and PTEN 0.543/0.550 ([same table](../results/metadata/task5_metrics.csv); sample unit: variants held out by protein). LOPO is limited to two Gate A PASS proteins and should not be generalized beyond this pair.

Paired B4−B2 AP effects are −0.087 [−0.163, −0.012] for CYP2C9 within-protein, −0.021 [−0.116, 0.071] for PTEN within-protein, −0.004 [−0.036, 0.025] for CYP2C9 held-out, and +0.007 [−0.041, 0.050] for PTEN held-out ([paired bootstrap summary](../results/metadata/task5_paired_bootstrap_summary.csv); sample units: variants and measured residue positions). The held-out intervals include zero. These are comparisons of B2/B4, which use ESM1v; they do not test whether ESM2 equals either learned feature set.

## Structure and sensitivity

Experimental structures mapped to 92.45% of CYP2C9 reference positions and 76.18% of PTEN reference positions, with identity 0.9978 and 1.0000 ([structure provenance](../results/metadata/task5_provenance.json); sample unit: reference residue positions). Missing residues remain missing. B3/B4 preprocessing and imputation are fit inside training folds ([structure metrics](../results/metadata/task5_metrics.csv); sample unit: variants).

The canonical sensitivity table retains controls rather than selecting only favorable results. Function-assay score-source checks are PASS with AP delta 0.000000 for both proteins/models; complete-case structure checks retain 637/694 CYP2C9 and 839/1,361 PTEN variants, with B4 AP deltas −0.004537 and +0.052843; label-confidence and assay-system comparisons are NOT_APPLICABLE because regimes/systems are confounded ([sensitivity results](../results/metadata/task5_sensitivity_results.csv); sample units: variants, positions, or proteins as named by each row). Proximity rows are descriptive and do not establish a mechanism. Figures 5–6 summarize these comparisons and deterministically selected post-analysis cases using a recorded rule.
