# Novelty Audit

The descriptive comparison of protein abundance and functional readouts is not new. This project therefore treats paired-assay harmonization as infrastructure and asks a narrower benchmarking question: do pretrained variant-effect scores recover functional impairment among variants whose measured abundance is preserved, and does that behavior generalize across proteins?

Cagiada et al. (2021) is the direct precedent for paired abundance/activity phenotype classes and structure-aware interpretation. Liao and Lehner (2026) is a more recent precedent for activity–abundance residuals, structure-informed analyses, and variant-predictor comparisons. This repository does **not** claim first discovery of abundance-preserved functional loss, allostery, or residual biology. Its intended contribution is a leakage-controlled, phenotype-conditioned benchmark with explicit score-source audits, position-cluster uncertainty, and a cross-protein extension; if that distinction is not sustained by the final evidence, the work should be labeled a reproducibility/extension study.

| Study | Proteins relevant here | Endpoint | Models | Validation split | Uncertainty and confounds | What this project adds |
|---|---|---|---|---|---|---|
| Cagiada et al. (2021) | PTEN, NUDT15 | Joint abundance/activity classes | Structural and sequence analyses | Descriptive | Compared different cellular assays | Not a first claim; supplies prior mechanism-oriented analysis |
| Amorosi et al. (2021) | CYP2C9 | Abundance and probe-binding activity | Variant predictors as supporting analyses | Within one protein | Human-cell abundance versus yeast activity | Adds a common cross-protein, phenotype-conditioned benchmark |
| Chiasson et al. (2020) | VKORC1 | Abundance and carboxylation activity | Supporting predictors | Within one protein | Both assays in human cells but with different reporters | Adds grouped validation and cross-protein comparison |
| Gersing et al. (2023/2024) | GCK | Yeast complementation and abundance-PCA | Stability/conservation analyses | Within one protein | Assays share yeast context but use different selection mechanisms | Adds a common score audit and held-out-protein evaluation |
| Liao and Lehner (2026) | PTEN, GCK and broader panels | Activity–abundance residual and allosteric effects | ESM-1v and ThermoMPNN | Primarily comparative/structural | Separates predicted stability from total fitness | Makes leakage control, phenotype-stratified AP, score coverage and position-cluster uncertainty the primary benchmark |

The intended contribution remains an extension or reproducibility benchmark unless the completed analysis demonstrates a distinction not already covered by these studies. Claims will be rewritten around the final evidence, including null or heterogeneous results.

## Sources

- ProteinGym reference metadata and distribution: <https://github.com/OATML-Markslab/ProteinGym>
- Cagiada et al.: <https://academic.oup.com/mbe/article/38/8/3235/6199445>
- Amorosi et al.: <https://pmc.ncbi.nlm.nih.gov/articles/PMC8456167/>
- Chiasson et al.: <https://elifesciences.org/articles/58026>
- Gersing activity map: <https://pmc.ncbi.nlm.nih.gov/articles/PMC10131484/>
- Gersing abundance map: <https://pmc.ncbi.nlm.nih.gov/articles/PMC11021015/>
- Liao and Lehner: <https://www.nature.com/articles/s41467-026-74517-8>
