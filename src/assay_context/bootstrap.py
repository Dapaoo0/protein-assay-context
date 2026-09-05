"""Position-clustered bootstrap for paired model comparisons."""

from __future__ import annotations

from collections.abc import Sequence

import numpy as np
import pandas as pd

from .evaluate import evaluate_binary


def _labels(frame: pd.DataFrame) -> np.ndarray:
    if "y" in frame:
        return pd.to_numeric(frame["y"], errors="raise").to_numpy(dtype=int)
    if "label" in frame:
        return pd.to_numeric(frame["label"], errors="raise").to_numpy(dtype=int)
    if "primary_endpoint" in frame:
        return frame["primary_endpoint"].eq("function-impaired").astype(int).to_numpy()
    raise ValueError("bootstrap frame needs y, label, or primary_endpoint")


def cluster_bootstrap(
    frame: pd.DataFrame,
    seed: int,
    n_boot: int,
    score_columns: Sequence[str] | None = None,
) -> pd.DataFrame:
    """Resample residue positions with replacement, retaining all their rows.

    A single sampled-position sequence is reused for every score column, so
    model contrasts are paired.  The row and position tuples are retained in
    output to make the cluster-preservation contract auditable.
    """
    if not isinstance(frame, pd.DataFrame):
        raise TypeError("frame must be a pandas DataFrame")
    if "position" not in frame:
        raise ValueError("bootstrap frame is missing required column: position")
    if not isinstance(n_boot, (int, np.integer)) or n_boot <= 0:
        raise ValueError("n_boot must be a positive integer")
    labels = _labels(frame)
    if score_columns is None:
        ignored = {"position", "variant", "y", "label", "primary_endpoint", "protein"}
        score_columns = [column for column in frame.columns if column not in ignored and pd.api.types.is_numeric_dtype(frame[column])]
    score_columns = list(score_columns)
    missing = sorted(set(score_columns) - set(frame.columns))
    if missing:
        raise ValueError(f"bootstrap frame is missing score columns: {', '.join(missing)}")
    if not score_columns:
        raise ValueError("at least one score column is required")
    positions = frame["position"].drop_duplicates().to_numpy()
    row_numbers = np.arange(len(frame), dtype=int)
    groups = [row_numbers[frame["position"].eq(position).to_numpy()] for position in positions]
    score_arrays = {model: pd.to_numeric(frame[model], errors="coerce").to_numpy(dtype=float) for model in score_columns}
    rng = np.random.default_rng(seed)
    rows: list[dict[str, object]] = []
    for draw in range(n_boot):
        selected = rng.integers(0, len(groups), size=len(groups))
        sampled_positions = tuple(positions[index].item() if hasattr(positions[index], "item") else positions[index] for index in selected)
        sampled_indices = tuple(int(value) for index in selected for value in groups[index])
        # Position resampling can contain a cluster more than once.  This is
        # represented explicitly by the repeated row indices in the tuple.
        sampled_y = labels[np.asarray(sampled_indices, dtype=int)]
        for model in score_columns:
            values = score_arrays[model][np.asarray(sampled_indices, dtype=int)]
            result = evaluate_binary(sampled_y, values) if np.isfinite(values).all() else {
                "n_variants": float(len(values)), "n_positives": float(sampled_y.sum()),
                "n_negatives": float(len(values) - sampled_y.sum()), "prevalence": float(sampled_y.mean()),
                "positives": float(sampled_y.sum()), "negatives": float(len(values) - sampled_y.sum()),
                "ap": np.nan, "ap_prevalence": np.nan, "ap_lift": np.nan, "auroc": np.nan, "valid": False,
            }
            rows.append({"draw": draw, "model": model, **result,
                         "sampled_positions": sampled_positions,
                         "sampled_row_indices": sampled_indices,
                         "invalid_reason": None if result["valid"] else ("single_class" if result["n_positives"] in (0, result["n_variants"]) else "non_finite_score"),
                         "seed": int(seed)})
    return pd.DataFrame(rows)


def paired_bootstrap_contrast(
    frame: pd.DataFrame,
    seed: int,
    n_boot: int,
    *,
    score_a: str,
    score_b: str,
) -> pd.DataFrame:
    """Compute paired AP contrast with invalid fraction over all draws.

    A draw is valid only when both model metrics are valid and non-missing.
    Invalid draws remain represented in the returned table; the reported
    fraction is computed before any CI filtering.
    """
    draws = cluster_bootstrap(frame, seed, n_boot, score_columns=[score_a, score_b])
    wide = draws.pivot(index="draw", columns="model", values=["ap", "valid", "sampled_positions", "sampled_row_indices"]).reset_index()
    wide.columns = ["_".join(str(part) for part in column if part) for column in wide.columns]
    ap_a, ap_b = f"ap_{score_a}", f"ap_{score_b}"
    valid_a, valid_b = f"valid_{score_a}", f"valid_{score_b}"
    wide["paired_ap_contrast"] = wide[ap_b] - wide[ap_a]
    wide["valid"] = wide[valid_a].astype(bool) & wide[valid_b].astype(bool) & wide[ap_a].notna() & wide[ap_b].notna()
    wide.loc[~wide["valid"], "paired_ap_contrast"] = np.nan
    invalid_fraction = float((~wide["valid"]).sum() / len(wide))
    wide["invalid_fraction"] = invalid_fraction
    wide["unstable"] = invalid_fraction > 0.10
    wide["seed"] = int(seed)
    wide["n_boot"] = int(n_boot)
    return wide


def position_block_permutation(
    frame: pd.DataFrame,
    seed: int,
    n_permutations: int,
    score_columns: Sequence[str],
) -> pd.DataFrame:
    """Observed-cohort null by shuffling complete label vectors by position.

    Positions are stratified by their number of measured variants.  Within each
    equal-size stratum, complete vectors of labels are permuted across score
    positions, preserving every vector's within-position heterogeneity and all
    group sizes.  This is a control for position-linked label structure, not a
    claim of independent variant-level exchangeability.
    """
    if "position" not in frame:
        raise ValueError("permutation frame is missing required column: position")
    labels = _labels(frame)
    score_columns = list(score_columns)
    missing = sorted(set(score_columns) - set(frame.columns))
    if missing:
        raise ValueError(f"permutation frame is missing score columns: {', '.join(missing)}")
    if n_permutations <= 0:
        raise ValueError("n_permutations must be positive")
    positions = frame["position"].drop_duplicates().tolist()
    groups = {position: np.flatnonzero(frame["position"].to_numpy() == position) for position in positions}
    strata: dict[int, list[object]] = {}
    for position, indices in groups.items():
        strata.setdefault(len(indices), []).append(position)
    rng = np.random.default_rng(seed)
    scores = {model: pd.to_numeric(frame[model], errors="coerce").to_numpy(dtype=float) for model in score_columns}
    rows: list[dict[str, object]] = []
    rule = "shuffle complete label vectors within equal-size residue-position strata"
    for permutation in range(n_permutations):
        shuffled = labels.copy()
        for indices_in_stratum in strata.values():
            source_vectors = [labels[groups[position]] for position in indices_in_stratum]
            for target_position, source_index in zip(indices_in_stratum, rng.permutation(len(indices_in_stratum)), strict=True):
                shuffled[groups[target_position]] = source_vectors[source_index]
        for model in score_columns:
            values = scores[model]
            result = evaluate_binary(shuffled, values) if np.isfinite(values).all() else {
                "n_variants": float(len(values)), "n_positives": float(shuffled.sum()), "n_negatives": float(len(values) - shuffled.sum()),
                "positives": float(shuffled.sum()), "negatives": float(len(values) - shuffled.sum()), "prevalence": float(shuffled.mean()),
                "ap": np.nan, "ap_prevalence": np.nan, "ap_lift": np.nan, "auroc": np.nan, "valid": False,
            }
            rows.append({"permutation": permutation, "model": model, **result, "seed": int(seed), "n_permutations": int(n_permutations), "label_block_rule": rule})
    return pd.DataFrame(rows)


observed_cohort_permutation = position_block_permutation
