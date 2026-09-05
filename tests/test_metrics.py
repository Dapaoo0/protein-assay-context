from __future__ import annotations

import importlib.util
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from assay_context.evaluate import (
    assert_common_cohort,
    evaluate_binary,
    evaluate_primary,
    substitution_descriptors,
)

_SCRIPT_SPEC = importlib.util.spec_from_file_location("benchmark", Path("scripts/04_benchmark.py"))
assert _SCRIPT_SPEC and _SCRIPT_SPEC.loader
_benchmark = importlib.util.module_from_spec(_SCRIPT_SPEC)
_SCRIPT_SPEC.loader.exec_module(_benchmark)


def test_perfect_ranking_has_unit_average_precision() -> None:
    result = evaluate_binary(np.array([0, 0, 1, 1]), np.array([0.1, 0.2, 0.8, 0.9]))
    assert result["ap"] == 1.0
    assert result["auroc"] == 1.0
    assert result["prevalence"] == 0.5


def test_single_class_metrics_are_flagged_invalid() -> None:
    result = evaluate_binary(np.array([1, 1, 1]), np.array([0.1, 0.2, 0.3]))
    assert result["valid"] is False
    assert np.isnan(result["ap"])


def test_primary_evaluation_uses_only_primary_eligible_rows_and_macro_is_two_proteins() -> None:
    frame = pd.DataFrame(
        {
            "protein": ["CYP2C9", "CYP2C9", "PTEN", "PTEN", "GCK"],
            "variant": ["A1V", "A1L", "A1V", "A1L", "A1V"],
            "position": [1, 1, 1, 1, 1],
            "wt": ["A"] * 5,
            "mutant": ["V", "L", "V", "L", "V"],
            "primary_eligible": [True, True, True, True, False],
            "primary_endpoint": ["function-preserved", "function-impaired", "function-impaired", "function-preserved", "not_eligible"],
            "ESM1v_single_impairment": [0.1, 0.9, 0.8, 0.2, 100.0],
            "ESM2_650M_impairment": [0.2, 0.8, 0.7, 0.3, 100.0],
        }
    )
    result = evaluate_primary(frame, bootstrap_draws=5)
    assert set(result.loc[result["protein"].isin(["CYP2C9", "PTEN"]), "n_variants"]) == {2.0}
    assert result.loc[result["protein"].eq("macro-average (2 proteins)"), "n_proteins"].eq(2).all()
    assert "GCK" not in set(result["protein"])


def test_common_cohort_rejects_different_canonical_ids() -> None:
    frame = pd.DataFrame(
        {
            "protein": ["P", "P"], "variant": ["A1V", "A1L"],
            "position": [1, 1], "wt": ["A", "A"], "mutant": ["V", "L"],
            "ESM1v_single_impairment": [0.1, 0.2], "ESM2_650M_impairment": [0.1, np.nan],
        }
    )
    with pytest.raises(ValueError, match="canonical IDs differ"):
        assert_common_cohort(frame)


def test_b1_descriptors_use_declared_named_scales() -> None:
    frame = pd.DataFrame({"wt": ["D", "T", "I", "K", "R", "H"], "mutant": ["K", "A", "K", "R", "K", "T"]})
    result = substitution_descriptors(frame)
    # Kyte-Doolittle hydropathy; acidic/basic charge convention is explicit.
    assert result["charge_change"].tolist() == [2, 0, 1, 0, 0, -1]
    assert np.allclose(result["hydrophobicity_change"], [-0.4, 2.5, -8.4, -0.6, 0.6, 2.5])


def test_primary_reader_rejects_same_token_with_different_frozen_key(tmp_path: Path) -> None:
    cohort = pd.DataFrame({
        "protein": ["CYP2C9"], "protein_id": ["CP2C9_HUMAN"], "sequence_hash": ["a" * 64],
        "variant": ["A1V"], "wt": ["A"], "position": [1], "mutant": ["V"],
        "primary_eligible": [True], "primary_endpoint": ["function-impaired"],
    })
    atlas = pd.DataFrame({
        "protein": ["CYP2C9"], "protein_id": ["CP2C9_HUMAN"], "sequence_hash": ["b" * 64],
        "variant": ["A1V"], "wt": ["A"], "position": [1], "mutant": ["V"],
        "ESM1v_single_impairment": [1.0], "ESM2_650M_impairment": [1.0],
    })
    cohort_path, atlas_path = tmp_path / "cohort.csv", tmp_path / "atlas.csv"
    cohort.to_csv(cohort_path, index=False)
    atlas.to_csv(atlas_path, index=False)
    with pytest.raises(ValueError, match="canonical|sequence_hash"):
        _benchmark._read_primary(cohort_path, atlas_path)


def test_frozen_input_hash_mismatch_is_rejected(tmp_path: Path) -> None:
    cohort_path, atlas_path, provenance_path = tmp_path / "cohort", tmp_path / "atlas", tmp_path / "provenance.json"
    cohort_path.write_bytes(b"cohort-v1")
    atlas_path.write_bytes(b"atlas-v1")
    provenance_path.write_text('{"input_sha256": {"cohort": "old", "atlas": "old"}}', encoding="utf-8")
    with pytest.raises(ValueError, match="input SHA256 mismatch"):
        _benchmark.validate_frozen_inputs(cohort_path, atlas_path, provenance_path)
