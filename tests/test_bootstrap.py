from __future__ import annotations

import numpy as np
import pandas as pd

from assay_context.bootstrap import (
    cluster_bootstrap,
    paired_bootstrap_contrast,
    position_block_permutation,
)


def _frame() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "variant": [f"A{position}{mutant}" for position in range(1, 7) for mutant in ("V", "L")],
            "position": [position for position in range(1, 7) for _ in ("V", "L")],
            "y": [0, 1] * 6,
            "model_a": np.linspace(0, 1, 12),
            "model_b": np.linspace(1, 0, 12),
        }
    )


def test_cluster_bootstrap_preserves_whole_positions_and_is_reproducible() -> None:
    frame = _frame()
    result = cluster_bootstrap(frame, seed=2026, n_boot=20)
    again = cluster_bootstrap(frame, seed=2026, n_boot=20)
    pd.testing.assert_frame_equal(result, again)
    for sampled in result["sampled_positions"]:
        assert set(sampled).issubset(set(frame["position"]))
    for sampled_rows in result["sampled_row_indices"]:
        for position in set(frame.iloc[list(sampled_rows)]["position"]):
            assert set(frame.index[frame["position"].eq(position)]).issubset(set(sampled_rows))


def test_cluster_bootstrap_flags_single_class_draws() -> None:
    frame = _frame().assign(y=1)
    result = cluster_bootstrap(frame, seed=2026, n_boot=4)
    assert result["valid"].eq(False).all()
    assert result["invalid_reason"].eq("single_class").all()


def test_paired_contrast_counts_mixed_validity_before_dropna() -> None:
    frame = _frame().copy()
    frame.loc[frame["position"].eq(3), "model_b"] = np.nan
    result = paired_bootstrap_contrast(frame, seed=2026, n_boot=20, score_a="model_a", score_b="model_b")
    assert result["draw"].nunique() == 20
    assert result["valid"].sum() < 20
    assert result["invalid_fraction"].iloc[0] > 0


def test_position_block_permutation_is_deterministic_and_preserves_label_vectors() -> None:
    frame = _frame()
    result = position_block_permutation(frame, seed=2026, n_permutations=10, score_columns=["model_a"])
    again = position_block_permutation(frame, seed=2026, n_permutations=10, score_columns=["model_a"])
    pd.testing.assert_frame_equal(result, again)
    assert result["permutation"].nunique() == 10
    assert result["n_variants"].eq(len(frame)).all()
    assert result["label_block_rule"].eq("shuffle complete label vectors within equal-size residue-position strata").all()
