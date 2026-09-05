from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from assay_context.sensitivity import (
    assert_common_variant_ids,
    assert_complete_b4_b2_contrasts,
    build_structure_features,
    compare_structure_models,
    paired_prediction_bootstrap,
    sensitivity_statuses,
    structure_complete_case,
    validate_sensitivity_evidence,
    verify_external_sha256,
)


def _frame() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "protein": ["CYP2C9"] * 8,
            "variant": [f"A{p}{m}" for p, m in zip((1, 1, 2, 2, 3, 3, 4, 4), ("V", "L") * 4, strict=True)],
            "position": [1, 1, 2, 2, 3, 3, 4, 4],
            "wt": ["A"] * 8,
            "mutant": ["V", "L"] * 4,
            "primary_endpoint": ["function-impaired", "function-preserved"] * 4,
            "primary_eligible": [True] * 8,
            "ESM1v_single_impairment": np.linspace(0.1, 0.8, 8),
            "accessibility": [0.1, np.nan, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8],
            "contacts": [2, 3, 4, 5, 6, 7, 8, 9],
            "secondary_structure": ["H", "H", "E", "E", "C", "C", "H", "H"],
            "site_distance": [4.0, 5.0, 6.0, 7.0, np.nan, 9.0, 10.0, 11.0],
        }
    )


def test_common_cohort_rejects_different_canonical_ids() -> None:
    frame = _frame()
    frame.loc[0, "variant"] = "A1W"
    with pytest.raises(ValueError, match="canonical|cohort"):
        assert_common_variant_ids(frame, frame.iloc[1:])


def test_structure_features_preserve_missing_indicators() -> None:
    features = build_structure_features(_frame())
    assert "accessibility_missing" in features
    assert bool(features.loc[0, "accessibility_missing"]) is False
    assert bool(features.loc[1, "accessibility_missing"]) is True


def test_structure_models_reuse_frozen_fold_column_and_return_b3_b4() -> None:
    frame = _frame()
    frame["task4_fold"] = [0, 0, 1, 1, 2, 2, 3, 3]
    result = compare_structure_models(frame, fold_column="task4_fold")
    assert set(result["feature_set"]) == {"B3", "B4"}
    assert result.groupby("feature_set").size().to_dict() == {"B3": 8, "B4": 8}
    assert result.groupby("feature_set")["fold"].apply(list).to_dict()["B3"] == [0, 0, 1, 1, 2, 2, 3, 3]


def test_sensitivity_table_never_omits_frozen_controls() -> None:
    statuses = sensitivity_statuses()
    required = {"function_score_source", "label_confidence", "coverage", "structure_missingness", "study_system", "proximity_4A", "proximity_5A", "proximity_8A"}
    assert required.issubset(set(statuses["control"]))
    assert statuses["status"].isin({"PASS", "NOT_APPLICABLE", "FAILED"}).all()


def test_external_source_hash_is_fail_closed(tmp_path) -> None:
    source = tmp_path / "source.csv"
    source.write_bytes(b"frozen")
    import hashlib

    expected = hashlib.sha256(b"frozen").hexdigest()
    assert verify_external_sha256(str(source), expected)
    assert not verify_external_sha256(str(source), "0" * 64)


def test_paired_prediction_bootstrap_uses_shared_positions() -> None:
    base = _frame().assign(y_true=_frame()["primary_endpoint"].eq("function-impaired").astype(int))
    frame = pd.concat([base.assign(feature_set=model, prediction=base["ESM1v_single_impairment"] + offset) for model, offset in (("B2", 0), ("B3", 0.1), ("B4", 0.2))], ignore_index=True)
    draws = paired_prediction_bootstrap(frame, setting="within", seed=2026, n_boot=5)
    assert set(draws["model"]) == {"B2", "B3", "B4"}
    assert draws["draw"].nunique() == 5
    assert draws.groupby("draw")["sampled_positions"].nunique().eq(1).all()


def test_sensitivity_pass_requires_numeric_evidence() -> None:
    bad = pd.DataFrame([{"control": "x", "status": "PASS", "reason": "claimed"}])
    with pytest.raises(ValueError, match="numeric evidence"):
        validate_sensitivity_evidence(bad)


def test_insufficient_n_is_explicit_nonpass_with_reason() -> None:
    result = pd.DataFrame([{
        "control": "proximity_4A", "status": "FAILED", "reason": "insufficient class count",
        "n_variants": 2, "n_positions": 1, "metric": np.nan, "effect": np.nan,
        "evidence": "D:/protein_assay_work/release_artifacts/task5_structure_mapping.csv",
    }])
    validate_sensitivity_evidence(result)
    assert result.loc[0, "status"] in {"FAILED", "NOT_APPLICABLE"}
    assert "insufficient" in result.loc[0, "reason"]


def test_figure5_contrast_metadata_has_exact_four_strata() -> None:
    summary = pd.DataFrame([
        {"setting": setting, "protein": protein, "model": "B4-B2"}
        for setting in ("within-protein", "held-out protein")
        for protein in ("CYP2C9", "PTEN")
    ])
    assert_complete_b4_b2_contrasts(summary)
    with pytest.raises(ValueError, match="required four"):
        assert_complete_b4_b2_contrasts(summary.iloc[:-1])


def test_complete_case_excludes_any_missing_structure_predictor() -> None:
    complete = structure_complete_case(_frame())
    assert len(complete) == 6
    assert complete["accessibility"].notna().all()
