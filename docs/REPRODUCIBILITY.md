# Reproducibility

## Offline environment and CI

Use Python 3.11 with the committed lockfile:

```powershell
uv sync --frozen
uv run python -m pytest -q
uv run python -m ruff check src tests scripts
uv build
```

The test suite is deterministic and offline after dependencies are installed. CI runs the same commands on Windows and Ubuntu with Python 3.11 ([workflow](../.github/workflows/ci.yml)). `tests/test_reproducibility.py` exercises a synthetic pair → label → canonical score join → grouped evaluation → position bootstrap twice with seed `2026`, and checks that cohort/metrics are byte-stable in memory.

## Full data pipeline

Raw inputs, archives, model scores, author tables, structures, and complete row-level outputs are stored on D: under `D:\protein_assay_work`. They are not required by the offline tests and are not redistributed. The commands below are the reproducible Windows invocation; use `/` in paths if running under a POSIX shell. Run the stages in order after downloading the official sources and checking every frozen SHA256.

### Frozen acquisition manifest

The source revisions and local destination filenames are fixed here. The author repositories are pinned to commits rather than mutable branches.

| Role | Frozen URL | Destination | SHA256 |
|---|---|---|---|
| ProteinGym reference metadata | [`144fe22`](https://raw.githubusercontent.com/OATML-Markslab/ProteinGym/144fe22b07dfaeec2b366f2346203a9838a55b4c/reference_files/DMS_substitutions.csv) | `D:\protein_assay_work\downloads\ProteinGym_DMS_substitutions.csv` | `a8f498011532a74aa9fe556a50555a75e928c5837d19c06a87592ae04049b308` |
| ProteinGym v1.3 archive | [ProteinGym v1.3](https://marks.hms.harvard.edu/proteingym/ProteinGym_v1.3/DMS_ProteinGym_substitutions.zip) | `D:\protein_assay_work\downloads\DMS_ProteinGym_substitutions_v1.3.zip` | `3a83766254ac9ac9984ec25cb73c6e010ea4418f5e35f143933e6b6e6473b921` |
| ProteinGym v1.3 zero-shot scores | [ProteinGym v1.3 zero-shot archive](https://marks.hms.harvard.edu/proteingym/ProteinGym_v1.3/zero_shot_substitutions_scores.zip) | `D:\protein_assay_work\downloads\zero_shot_substitutions_scores.zip` | `3fd7cdb5e78f1d43cabfabfeb6578c252b63af23ba2ab44db0094dc3a42de36d` |
| CYP2C9 author source (`data/CYP2C9_activity_abundance_scores.csv`, `e727ffe`) | [`e727ffe`](https://raw.githubusercontent.com/dunhamlab/CYP2C9/e727ffe1e4746311d44c09b4f6ffa35ca3acf665/data/CYP2C9_activity_abundance_scores.csv) | `D:\protein_assay_work\author_source\CYP2C9_activity_abundance_scores.csv` | `283883889947e47c47b08f2f096feaf7ad5dd49d0f79432b8675c10150b61124` |
| PTEN abundance (`output_datatables/Composite_abundance_data.tsv`, `73443cf`) | [`73443cf`](https://raw.githubusercontent.com/matreyeklab/pten_composite/73443cfa612ae6e570347f96461ed71a2641608a/output_datatables/Composite_abundance_data.tsv) | `D:\protein_assay_work\author_source\PTEN_Composite_abundance_data.tsv` | `e8172400de12c53043556f3d0b07e44e5a23c21582b14f3bc8f6d3dad9c0fe6c` |
| PTEN function (`input_datatables/Mighell_PTEN_phosphatase_data.csv`, `73443cf`) | [`73443cf`](https://raw.githubusercontent.com/matreyeklab/pten_composite/73443cfa612ae6e570347f96461ed71a2641608a/input_datatables/Mighell_PTEN_phosphatase_data.csv) | `D:\protein_assay_work\author_source\PTEN_Mighell_phosphatase.csv` | `6d0bc5e0f6a50f59ca272a19b5756107a089473d2d8e54a706e204da27601ecd` |
| GCK abundance (`data/abundance.csv`, `e77eeb3`) | [`e77eeb3`](https://raw.githubusercontent.com/KULL-Centre/_2024_Gersing_GCKabundance/e77eeb3a3a3ec89baf4ac0e4196a041f7d5d32ae/data/abundance.csv) | `D:\protein_assay_work\author_source\GCK_abundance.csv` | `75fdeece9cc1c1a35a417956b1e9d5c62c39fe0a4a7c91a364bc43111c837dc7` |
| GCK function (`data/activity.csv`, `e77eeb3`) | [`e77eeb3`](https://raw.githubusercontent.com/KULL-Centre/_2024_Gersing_GCKabundance/e77eeb3a3a3ec89baf4ac0e4196a041f7d5d32ae/data/activity.csv) | `D:\protein_assay_work\author_source\GCK_activity.csv` | `f7efb1d51dd380fcd275c56da70592634b945579bbdf3acc065ffe226c19e7cd` |
| VKORC1 combined author source (`f80bb91`) | [`f80bb91`](https://raw.githubusercontent.com/FowlerLab/VKOR/f80bb91695967393adacd7b33ded6632283d42f1/vkor_abundance_activity_scores_variants.csv) | `D:\protein_assay_work\author_source\VKOR_combined.csv` | `ea4bf42a40f1a7965667eca37318b2f17fc34ae8174f656939d5626e56523a8b` |
| CYP2C9 structure | [`1R9O.pdb`](https://files.rcsb.org/download/1R9O.pdb) | `D:\protein_assay_work\structures\1R9O.pdb` | `16860aca35da55e576db41e3ff6a69227b41f771ba5ff4691d1d61eeef10e136` |
| PTEN structure | [`1D5R.pdb`](https://files.rcsb.org/download/1D5R.pdb) | `D:\protein_assay_work\structures\1D5R.pdb` | `efd284b66f0f66d6c78c42aebe743c2095165878c097eff1eddee1f57653914c` |
| CYP2C9 UniProt site annotation | [`P11712.txt`](https://rest.uniprot.org/uniprotkb/P11712.txt) | `D:\protein_assay_work\structures\P11712.txt` | `04225b238d8651585744c51ef485bf68776794b5c950c8496e6597f12766b9cb` |
| PTEN UniProt site annotation | [`P60484.txt`](https://rest.uniprot.org/uniprotkb/P60484.txt) | `D:\protein_assay_work\structures\P60484.txt` | `e0f06e38d1d104d40bba113d5d1d6d904c5f51d68ff4fdca2701f263505d4c75` |

For each download, use the exact destination and verify before parsing. This example is executable in PowerShell (the URL list above is authoritative):

```powershell
$work = 'D:\protein_assay_work'
$url = 'https://raw.githubusercontent.com/OATML-Markslab/ProteinGym/144fe22b07dfaeec2b366f2346203a9838a55b4c/reference_files/DMS_substitutions.csv'
$dest = Join-Path $work 'downloads\ProteinGym_DMS_substitutions.csv'
Invoke-WebRequest -Uri $url -OutFile $dest
Get-FileHash -Algorithm SHA256 -LiteralPath $dest
```

Download the v1.3 archive to the manifest destination and expand it to `D:\protein_assay_work\assays`; retain the archive for the checksum audit. Do not commit or redistribute the archive, raw author files, PDBs, or extracted row-level data.

### Exact stage commands

```powershell
uv run python scripts/01_inventory.py `
  --reference D:/protein_assay_work/downloads/ProteinGym_DMS_substitutions.csv `
  --archive D:/protein_assay_work/downloads/DMS_ProteinGym_substitutions_v1.3.zip `
  --assay-dir D:/protein_assay_work/assays `
  --config config/candidates.yaml `
  --output-dir results/metadata

uv run python scripts/02_pair_assays.py `
  --config config/candidates.yaml `
  --inventory results/metadata/assay_inventory.csv `
  --assay-dir D:/protein_assay_work/assays `
  --author-dir D:/protein_assay_work/author_source `
  --processed-dir D:/protein_assay_work/processed `
  --output-dir results/metadata

uv run python scripts/03_build_atlas.py `
  --cohort D:/protein_assay_work/processed/paired_phenotype_cohort.csv `
  --score-dir D:/protein_assay_work/processed/model_scores `
  --archive D:/protein_assay_work/downloads/zero_shot_substitutions_scores.zip `
  --metadata-dir results/metadata `
  --figure-dir results/figures `
  --config config/analysis_freeze.yaml

uv run python scripts/04_benchmark.py `
  --cohort D:/protein_assay_work/processed/paired_phenotype_cohort.csv `
  --atlas D:/protein_assay_work/release_artifacts/assay_atlas.csv `
  --metadata-dir results/metadata `
  --figure-dir results/figures `
  --seed 2026 `
  --bootstrap-draws 2000 `
  --permutation-draws 2000

uv run python scripts/05_structure_sensitivity.py `
  --structure-dir D:/protein_assay_work/structures `
  --metadata-dir results/metadata `
  --figure-dir results/figures `
  --config config/analysis_freeze.yaml
```

The score archive is the frozen ProteinGym zero-shot archive listed above; its extracted model scores are local pipeline inputs. The acquisition manifest and `config/analysis_freeze.yaml` record the remaining revisions and hashes. If it is absent, download that exact URL and verify its hash before running stage 03.

`config/cohort_freeze.json` is the release contract: v0.1 contains the amended primary binary benchmark only. Spearman score/function and grouped-OOF isotonic-residual endpoints are explicitly deferred, so no secondary inference is claimed.

The atlas and row-level prediction/mapping outputs are written to `D:\protein_assay_work\release_artifacts` and are ignored by the public tree. Audit samples and case tables are also generated locally and ignored. Compact metrics, summaries, figures, and provenance remain tracked. The reproducibility test validates the release-artifact manifest hashes when that local directory is present.

The plan budget was working data ≤2 GB, peak environment/archive storage ≤8 GB, and RAM ≤8 GB; these are engineering targets, not claims that a fresh download will fit automatically. Cleanup is deferred and must follow a separate inventory/verification step.
