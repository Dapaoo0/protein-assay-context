# Data provenance and redistribution

## Frozen sources

| Role | Source and revision | Retrieval/hash evidence |
|---|---|---|
| Assay metadata and score distribution | [ProteinGym](https://github.com/OATML-Markslab/ProteinGym), revision `144fe22b07dfaeec2b366f2346203a9838a55b4c`; [commit-pinned reference CSV](https://raw.githubusercontent.com/OATML-Markslab/ProteinGym/144fe22b07dfaeec2b366f2346203a9838a55b4c/reference_files/DMS_substitutions.csv) | Metadata SHA256 `a8f498011532a74aa9fe556a50555a75e928c5837d19c06a87592ae04049b308`; score archive SHA256 `3fd7cdb5e78f1d43cabfabfeb6578c252b63af23ba2ab44db0094dc3a42de36d`, 1,911,045,703 bytes ([freeze](../config/analysis_freeze.yaml)) |
| CYP2C9 author abundance/function | [dunhamlab/CYP2C9](https://github.com/dunhamlab/CYP2C9), revision `e727ffe1e4746311d44c09b4f6ffa35ca3acf665` | `283883889947e47c47b08f2f096feaf7ad5dd49d0f79432b8675c10150b61124` ([Task 2 provenance](../results/metadata/task2_provenance.json)) |
| PTEN author abundance/function | [matreyeklab/pten_composite](https://github.com/matreyeklab/pten_composite), revision `73443cfa612ae6e570347f96461ed71a2641608a` | `e8172400de12c53043556f3d0b07e44e5a23c21582b14f3bc8f6d3dad9c0fe6c` and `6d0bc5e0f6a50f59ca272a19b5756107a089473d2d8e54a706e204da27601ecd` ([Task 2 provenance](../results/metadata/task2_provenance.json)) |
| GCK author abundance/function | [KULL-Centre/_2024_Gersing_GCKabundance](https://github.com/KULL-Centre/_2024_Gersing_GCKabundance), revision `e77eeb3a3a3ec89baf4ac0e4196a041f7d5d32ae` | `75fdeece9cc1c1a35a417956b1e9d5c62c39fe0a4a7c91a364bc43111c837dc7` and `f7efb1d51dd380fcd275c56da70592634b945579bbdf3acc065ffe226c19e7cd` ([Task 2 provenance](../results/metadata/task2_provenance.json)) |
| VKORC1 author abundance/function | [FowlerLab/VKOR](https://github.com/FowlerLab/VKOR), revision `f80bb91695967393adacd7b33ded6632283d42f1` | `ea4bf42a40f1a7965667eca37318b2f17fc34ae8174f656939d5626e56523a8b` ([Task 2 provenance](../results/metadata/task2_provenance.json)) |
| Experimental structures and independent sites | [RCSB 1R9O](https://files.rcsb.org/download/1R9O.pdb), [RCSB 1D5R](https://files.rcsb.org/download/1D5R.pdb), [UniProt P11712](https://rest.uniprot.org/uniprotkb/P11712.txt), [UniProt P60484](https://rest.uniprot.org/uniprotkb/P60484.txt) | Structure/site hashes, revisions, and mapping rules ([Task 5 provenance](../results/metadata/task5_provenance.json)) |

The GCK abundance hash above is recorded from the frozen manifest; the exact complete value is also machine-checked in `config/analysis_freeze.yaml`.

## Storage and release policy

Raw author/ProteinGym archives, raw PDB/UniProt files, model archives, the paired cohort, audit samples, case tables, and complete draw-level files remain local under `D:\protein_assay_work`. Seven complete row-level source-derived outputs are retained intact in `D:\protein_assay_work\release_artifacts` with SHA256, byte count, and row count in its local `manifest.json`. Compact aggregates, figures, and provenance remain tracked; no variant-level CSV is distributed in the public tree.

The ProteinGym code repository is [MIT](https://raw.githubusercontent.com/OATML-Markslab/ProteinGym/main/LICENSE). The GCK source repository is MIT; the VKOR source repository is BSD-3-Clause. No LICENSE file was found at the frozen CYP2C9 or PTEN source-repository paths during the release audit. A repository license does not unambiguously relicense third-party DMS data or derived archive contents. Accordingly, this project does not claim redistribution permission for raw or row-level source-derived data; the public release contains code and compact summaries only. The project code is MIT-licensed, while upstream data retain their original terms.
