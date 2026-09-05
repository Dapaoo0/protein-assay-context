"""Offline, deterministic release-contract checks."""

from __future__ import annotations

import hashlib
import json
import subprocess
import tarfile
from pathlib import Path

import numpy as np
import pandas as pd
from PIL import Image

from assay_context.artifacts import external_artifact_path
from assay_context.atlas import summarize_phenotypes
from assay_context.bootstrap import cluster_bootstrap
from assay_context.evaluate import evaluate_primary
from assay_context.scores import join_scores


def _synthetic_pairs() -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for position in range(1, 7):
        for mutant, impaired in (("V", position % 2 == 0), ("L", position % 2 == 1)):
            rows.append(
                {
                    "protein": "SYNTH",
                    "protein_id": "P00000",
                    "sequence_hash": "a" * 64,
                    "wt": "A",
                    "position": position,
                    "mutant": mutant,
                    "variant": f"A{position}{mutant}",
                    "abundance_label": "preserved",
                    "function_label": "impaired" if impaired else "preserved",
                    "primary_eligible": True,
                    "primary_endpoint": "function-impaired" if impaired else "function-preserved",
                }
            )
    return pd.DataFrame(rows)


def _run_synthetic() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    pairs = _synthetic_pairs()
    scores = pairs[["protein", "protein_id", "sequence_hash", "wt", "position", "mutant", "variant"]].copy()
    scores["ESM1v_single"] = np.linspace(0.01, 0.99, len(scores))
    scores["ESM2_650M"] = scores["ESM1v_single"] ** 1.1
    joined = join_scores(pairs, scores)
    joined["ESM1v_single_impairment"] = -joined["ESM1v_single"]
    joined["ESM2_650M_impairment"] = -joined["ESM2_650M"]
    direct = evaluate_primary(joined, seed=2026, bootstrap_draws=200)
    atlas = summarize_phenotypes(pairs, seed=2026, bootstrap_draws=200)
    draws = cluster_bootstrap(
        joined.assign(y=joined["primary_endpoint"].eq("function-impaired").astype(int)),
        seed=2026,
        n_boot=200,
        score_columns=["ESM1v_single_impairment", "ESM2_650M_impairment"],
    )
    return direct, atlas, draws


def test_synthetic_end_to_end_is_deterministic_and_position_clustered() -> None:
    first = _run_synthetic()
    second = _run_synthetic()
    for left, right in zip(first, second, strict=True):
        pd.testing.assert_frame_equal(left, right)
    assert first[0].loc[first[0]["protein"].eq("SYNTH"), "n_variants"].eq(12).all()
    assert first[1]["denominator"].eq(12).all()
    assert first[2]["draw"].nunique() == 200
    for sampled in first[2]["sampled_row_indices"]:
        assert len(sampled) == 12
        assert all(
            sampled[offset] // 2 == sampled[offset + 1] // 2
            for offset in range(0, len(sampled), 2)
        )


def test_tracked_compact_metrics_have_defined_denominators() -> None:
    root = Path(__file__).resolve().parents[1]
    for relative in ("results/metadata/task4_primary_metrics.csv", "results/metadata/task5_metrics.csv"):
        frame = pd.read_csv(root / relative)
        assert frame["n_variants"].dropna().gt(0).all(), relative
        assert frame["prevalence"].dropna().between(0, 1).all(), relative
        assert frame["ap"].dropna().between(0, 1).all(), relative


def test_manifest_hashes_externalized_row_level_outputs() -> None:
    manifest = Path("D:/protein_assay_work/release_artifacts/manifest.json")
    if not manifest.exists():
        return
    payload = json.loads(manifest.read_text(encoding="utf-8-sig"))
    for item in payload["files"]:
        path = manifest.parent / item["filename"]
        digest = hashlib.sha256(path.read_bytes()).hexdigest()
        assert digest == item["sha256"]


def test_figure_manifest_matches_pngs() -> None:
    root = Path(__file__).resolve().parents[1]
    manifest = pd.read_csv(root / "results/metadata/figure_manifest.csv")
    assert len(manifest) == 6
    for row in manifest.itertuples(index=False):
        path = root / row.path
        with Image.open(path) as image:
            assert (image.width, image.height) == (row.width_px, row.height_px)
        assert hashlib.sha256(path.read_bytes()).hexdigest() == row.sha256


def test_tracked_tree_excludes_local_data_and_secret_patterns() -> None:
    root = Path(__file__).resolve().parents[1]
    tracked = subprocess.check_output(["git", "ls-files"], cwd=root, text=True).splitlines()
    variant_level_tables = {
        "results/metadata/assay_atlas.csv",
        "results/metadata/task2_audit_sample.csv",
        "results/metadata/task4_lopo_predictions.csv",
        "results/metadata/task4_oof_predictions.csv",
        "results/metadata/task5_case_studies.csv",
        "results/metadata/task5_complete_case_oof_predictions.csv",
        "results/metadata/task5_lopo_predictions.csv",
        "results/metadata/task5_oof_predictions.csv",
        "results/metadata/task5_structure_mapping.csv",
    }
    assert variant_level_tables.isdisjoint(tracked)
    forbidden_suffixes = (".env", ".zip", ".tar.gz", ".pdb", ".cif")
    assert not any(
        path.startswith((".venv/", "data/raw/", "data/cache/")) or path.lower().endswith(forbidden_suffixes)
        for path in tracked
    )
    secret_pattern = __import__("re").compile(r"(?i)(api[_-]?key|secret|password)\s*[:=]\s*['\"]")
    for relative in tracked:
        if relative.endswith((".md", ".py", ".yml", ".yaml", ".toml", ".cff")):
            assert not secret_pattern.search((root / relative).read_text(encoding="utf-8")), relative


def test_row_level_outputs_resolve_outside_public_tree() -> None:
    assert external_artifact_path("assay_atlas.csv").as_posix() == "D:/protein_assay_work/release_artifacts/assay_atlas.csv"
    assert external_artifact_path("task4_oof_predictions.csv").as_posix() == "D:/protein_assay_work/release_artifacts/task4_oof_predictions.csv"


def test_sdist_excludes_internal_review_and_row_level_payloads(tmp_path: Path) -> None:
    root = Path(__file__).resolve().parents[1]
    output = tmp_path / "dist"
    subprocess.run(
        ["uv", "build", "--offline", "--sdist", "--out-dir", str(output)],
        cwd=root,
        check=True,
        capture_output=True,
        text=True,
    )
    archive = next(output.glob("*.tar.gz"))
    with tarfile.open(archive, "r:gz") as bundle:
        relative_names = [name.split("/", maxsplit=1)[1] for name in bundle.getnames()]

    forbidden_basenames = {
        "assay_atlas.csv",
        "task2_audit_sample.csv",
        "task4_lopo_predictions.csv",
        "task4_oof_predictions.csv",
        "task5_case_studies.csv",
        "task5_complete_case_oof_predictions.csv",
        "task5_lopo_predictions.csv",
        "task5_oof_predictions.csv",
        "task5_structure_mapping.csv",
    }
    assert not any(name.startswith(".superpowers/") for name in relative_names)
    assert forbidden_basenames.isdisjoint(Path(name).name for name in relative_names)


def test_cohort_freeze_records_scope_and_dated_secondary_amendment() -> None:
    root = Path(__file__).resolve().parents[1]
    freeze = json.loads((root / "config/cohort_freeze.json").read_text(encoding="utf-8"))
    assert freeze["release_version"] == "v0.1"
    assert freeze["primary_proteins"] == ["CYP2C9", "PTEN"]
    assert freeze["primary_cohort"]["n_variants"] == 2055
    assert freeze["primary_cohort"]["n_measured_residue_positions"] == 550
    amendment = freeze["post_analysis_amendment"]
    assert amendment["date"] == "2026-09-05"
    assert set(amendment["deferred_endpoints"]) == {
        "Spearman score/function",
        "grouped-OOF isotonic residual",
    }
    assert amendment["status"] == "deferred"


def test_case_docs_and_config_use_post_analysis_selection_wording() -> None:
    root = Path(__file__).resolve().parents[1]
    case_paths = [
        root / "README.md",
        root / "docs/RESULTS.md",
        root / "results/metadata/task5_provenance.json",
        root / "scripts/05_structure_sensitivity.py",
    ]
    case_text = "\n".join(path.read_text(encoding="utf-8") for path in case_paths)
    assert "predeclared" not in case_text.lower()
    freeze_text = (root / "config/analysis_freeze.yaml").read_text(encoding="utf-8")
    assert "case_rule: deterministically selected post-analysis using a recorded rule" in freeze_text
    assert "deterministically selected post-analysis using a recorded rule" in case_text
