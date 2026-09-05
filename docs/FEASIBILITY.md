# Feasibility Audit

Audit date: 2026-09-05. ProteinGym reference revision: `144fe22b07dfaeec2b366f2346203a9838a55b4c`. DMS archive: ProteinGym v1.3.

The mapping gate passes for all four candidate pairs. Counts below are exact intersections of canonical single amino-acid substitutions after checking both assays against identical ProteinGym target sequences. They are not the per-assay totals from the metadata table.

| Protein | Shared missense variants | Shared positions | Reference sequences identical | WT mismatches | Mapping gate |
|---|---:|---:|---|---:|---|
| CYP2C9 | 4,421 | 464 | yes | 0 | PASS |
| PTEN | 4,839 | 379 | yes | 0 | PASS |
| GCK | 8,255 | 463 | yes | 0 | PASS |
| VKORC1 | 692 | 151 | yes | 0 | PASS |

This does not yet pass the complete Gate A. The primary binary endpoint also requires at least 50 function-impaired and 50 function-preserved abundance-preserved variants, each spanning at least 10 positions. Those labels require author-calibrated thresholds and uncertainty data and are frozen in Task 2. Gate B for pretrained-score coverage is deferred to Task 3.

## Preliminary gate verdict

The machine-readable `results/metadata/feasibility_summary.json` uses only `PASS` or `REDUCE_SCOPE` statuses and includes a reason for every gate. Mapping is `PASS` for the four pairs above. Class balance is `REDUCE_SCOPE` until Task 2 freezes author-calibrated thresholds, uncertainty handling, and minimum class/position counts. Model-score coverage is `REDUCE_SCOPE` until Task 3 audits score coverage and freezes the score cohort. Therefore the preliminary overall verdict is `REDUCE_SCOPE`, with these reasons preserved in `overall_reasons` and the nested `gates` object.

## Reproducible invocation and ignored input layout

After placing the downloaded inputs in the ignored `data/raw/` tree, run:

```text
uv run python scripts/01_inventory.py \
  --reference data/raw/ProteinGym_DMS_substitutions.csv \
  --archive data/raw/ProteinGym_v1.3/DMS_ProteinGym_substitutions.zip \
  --assay-dir data/raw/ProteinGym_v1.3/assays \
  --config config/candidates.yaml \
  --output-dir results/metadata
```

Expected local-only input layout (all `data/raw/` content is ignored):

```text
data/raw/
├── ProteinGym_DMS_substitutions.csv
└── ProteinGym_v1.3/
    ├── DMS_ProteinGym_substitutions.zip
    └── assays/
        ├── <configured-abundance-assay>.csv
        └── <configured-function-assay>.csv
```

The script checks both download hashes, validates and resolves every pair before reading assay files, and writes the inventory, overlap table, and gate summary to `results/metadata/`.

## Assay interpretation

- CYP2C9 pairs human-cell VAMP-seq abundance with an activity-based probe measured by click-seq in engineered yeast. The functional readout includes probe binding and catalytic context; it is not a pure kinetic measurement.
- PTEN pairs abundance in a human-derived cell line with rescue of PI3K-induced toxicity in yeast, a growth proxy for lipid phosphatase activity.
- GCK pairs yeast abundance-PCA with glucose-growth complementation in a triple hexokinase-deletion yeast strain. Both are growth selections but report different biological quantities.
- VKORC1 pairs VAMP-seq with a cell-surface gamma-carboxylation reporter in engineered human cells.

Differences between assay systems are recorded as possible confounders. A discrepant score is not by itself evidence of allostery or clinical pathogenicity.

## Storage check

The v1.3 DMS archive is 43,021,128 bytes compressed. Only eight assay CSVs were extracted into the ignored work directory on drive D. The zero-shot score archive was not downloaded; its server-reported size is approximately 1.91 GB and will be considered only after the phenotype cohort is frozen.

## Task 2 phenotype freeze

The paired phenotype cohort was frozen on 2026-09-05 with fixed audit seed 2026.
Author labels and thresholds were applied before any model-score evaluation; the
ProteinGym median bins were not used as endpoint labels. Non-missense rows are
excluded from the paired missense cohort and counted in the attrition audit.

| Protein | Paired missense | Abundance-preserved | Function-impaired | Function-preserved | Impaired positions | Preserved positions | Gate A |
|---|---:|---:|---:|---:|---:|---:|---|
| CYP2C9 | 4,421 | 1,244 | 412 | 282 | 191 | 145 | PASS |
| PTEN | 4,839 | 1,545 | 212 | 1,149 | 101 | 270 | PASS |
| GCK | 8,255 | 1,871 | 460 | 45 | 149 | 41 | REDUCE_SCOPE |
| VKORC1 | 692 | 232 | 4 | 178 | 4 | 90 | REDUCE_SCOPE |

Gate A therefore supports the primary binary endpoint for CYP2C9 and PTEN.
GCK and VKORC1 remain descriptive/continuous for this endpoint because one
class or its position minimum is not met. Full row-level outputs are kept in
the ignored `D:\protein_assay_work\processed` directory; tracked summaries,
source hashes, and the 20-variant-per-protein audit are in
`results/metadata/task2_*.csv` and `task2_provenance.json`.
