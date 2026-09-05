"""Validation and left joins for ProteinGym zero-shot model scores."""

from __future__ import annotations

from collections.abc import Iterable

import pandas as pd

CANONICAL_KEYS = ("protein_id", "sequence_hash", "wt", "position", "mutant")
PAIR_KEYS = ("protein", "variant")
CONTEXT_KEYS = ("construct", "offset")


def _require(frame: pd.DataFrame, columns: Iterable[str], label: str) -> None:
    missing = sorted(set(columns) - set(frame.columns))
    if missing:
        raise ValueError(f"{label} is missing required columns: {', '.join(missing)}")


def _canonical_keys(frame: pd.DataFrame, label: str) -> pd.DataFrame:
    result = frame.copy()
    if "variant" not in result and {"wt", "position", "mutant"}.issubset(result.columns):
        result["variant"] = (
            result["wt"].astype("string")
            + result["position"].astype("string")
            + result["mutant"].astype("string")
        )
    _require(result, (*CANONICAL_KEYS, *PAIR_KEYS), label)
    result["protein"] = result["protein"].astype("string")
    result["variant"] = result["variant"].astype("string").str.strip().str.upper()
    result["protein_id"] = result["protein_id"].astype("string").str.strip()
    result["sequence_hash"] = result["sequence_hash"].astype("string").str.strip().str.lower()
    result["wt"] = result["wt"].astype("string").str.strip().str.upper()
    result["mutant"] = result["mutant"].astype("string").str.strip().str.upper()
    result["position"] = pd.to_numeric(result["position"], errors="raise").astype("int64")
    if "construct" in result:
        result["construct"] = result["construct"].astype("string").str.strip()
    if "offset" in result:
        result["offset"] = pd.to_numeric(result["offset"], errors="raise").astype("int64")
    if result[list(CANONICAL_KEYS) + ["protein", "variant"]].isna().any().any() or result["variant"].eq("").any():
        raise ValueError(f"{label} contains missing canonical keys")
    expected = result["wt"] + result["position"].astype("string") + result["mutant"]
    if not result["variant"].eq(expected).all():
        raise ValueError(f"{label} has inconsistent variant and parsed canonical components")
    return result


def _to_long(scores: pd.DataFrame) -> pd.DataFrame:
    """Accept the compact long contract or a ProteinGym wide score table."""
    frame = _canonical_keys(scores, "scores")
    if {"model", "score"}.issubset(frame.columns):
        result = frame.copy()
        result["score"] = pd.to_numeric(result["score"], errors="coerce")
        return result
    identifiers = {
        "protein", "protein_id", "sequence_hash", "variant", "wt", "position", "mutant", "assay_id",
        "construct", "offset", "source", "source_archive", "direction", "higher_is_impaired",
    }
    models = [column for column in frame.columns if column not in identifiers]
    if not models:
        raise ValueError("scores must contain model/score columns or at least one model column")
    result = frame.melt(
        id_vars=[column for column in frame.columns if column in identifiers],
        value_vars=models,
        var_name="model",
        value_name="score",
    )
    result["score"] = pd.to_numeric(result["score"], errors="coerce")
    return result


def _score_keys(context_keys: tuple[str, ...] = ()) -> list[str]:
    return [*CANONICAL_KEYS, *PAIR_KEYS, *context_keys, "model"]


def _require_unique_source(long: pd.DataFrame, context_keys: tuple[str, ...], label: str) -> None:
    if long.duplicated(_score_keys(context_keys)).any():
        raise ValueError(f"{label} contains duplicate canonical score keys")


def join_scores(pairs: pd.DataFrame, scores: pd.DataFrame) -> pd.DataFrame:
    """Left-join validated model scores without dropping score-missing pairs.

    Scores may be supplied long (``model``, ``score``) or in the native
    ProteinGym wide format. Every canonical model key must be unique; this
    function is intentionally a single-source production join.
    """
    pair_frame = _canonical_keys(pairs, "pairs")
    if pair_frame.duplicated(list(CANONICAL_KEYS)).any():
        raise ValueError("pairs contains duplicate canonical keys")
    long = _to_long(scores)
    context_keys = tuple(key for key in CONTEXT_KEYS if key in pair_frame or key in long)
    for key in context_keys:
        if key not in pair_frame or key not in long:
            raise ValueError(f"{key} mismatch: both pairs and scores must declare it")
    if context_keys:
        base_keys = [*CANONICAL_KEYS, *PAIR_KEYS]
        context_check = long.merge(
            pair_frame[base_keys + list(context_keys)],
            on=base_keys,
            how="inner",
            suffixes=("_score", "_pair"),
        )
        for key in context_keys:
            if not context_check[f"{key}_score"].eq(context_check[f"{key}_pair"]).all():
                raise ValueError(f"{key} mismatch between pairs and scores")
    pair_hashes = pair_frame.groupby("protein_id")["sequence_hash"].unique().to_dict()
    for protein_id, group in long.groupby("protein_id", sort=False):
        allowed = set(pair_hashes.get(protein_id, []))
        if allowed and not set(group["sequence_hash"]).issubset(allowed):
            raise ValueError(f"scores sequence_hash mismatch for protein_id {protein_id}")
    join_keys = [*CANONICAL_KEYS, *PAIR_KEYS, *context_keys]
    _require_unique_source(long, context_keys, "scores")
    wide = long.pivot(index=join_keys, columns="model", values="score").reset_index()
    wide.columns.name = None
    result = pair_frame.merge(wide, on=join_keys, how="left", validate="one_to_one")
    coverage: dict[str, dict[str, int | float]] = {}
    total = len(pair_frame)
    for model in long["model"].drop_duplicates().tolist():
        matched = int(result[model].notna().sum())
        coverage[str(model)] = {"matched": matched, "total": total, "fraction": matched / total if total else 0.0}
    result.attrs["coverage"] = coverage
    result.attrs["consistency"] = {
        "duplicate_keys_rejected": True,
        "groups_checked": len(long),
        "groups_checked_by_model": long.groupby("model").size().to_dict(),
    }
    result.attrs["models"] = long["model"].drop_duplicates().tolist()
    return result


def audit_cross_assay_scores(
    abundance: pd.DataFrame, function: pd.DataFrame
) -> tuple[dict[str, int], pd.DataFrame]:
    """Audit exact numeric agreement without using function scores in production.

    The two inputs must each contain unique canonical keys.  Every matched
    model/variant is classified as exact or non-identical; no tolerance or
    coalescing is applied.
    """
    abundance_long = _to_long(abundance)
    function_long = _to_long(function)
    if "assay_id" not in abundance_long or "assay_id" not in function_long:
        raise ValueError("cross-assay audit requires assay_id in both sources")
    context_keys = tuple(
        key for key in CONTEXT_KEYS if key in abundance_long or key in function_long
    )
    for key in context_keys:
        if key not in abundance_long or key not in function_long:
            raise ValueError(f"{key} mismatch: both audit sources must declare it")
    _require_unique_source(abundance_long, context_keys, "abundance source")
    _require_unique_source(function_long, context_keys, "function source")
    keys = _score_keys(context_keys)
    merged = abundance_long.merge(
        function_long,
        on=keys,
        how="outer",
        suffixes=("_abundance", "_function"),
        indicator=True,
    )
    matched = merged.loc[merged["_merge"].eq("both")].copy()
    both_missing = matched["score_abundance"].isna() & matched["score_function"].isna()
    exact = matched["score_abundance"].eq(matched["score_function"]) | both_missing
    nonidentical = matched.loc[~exact].copy()
    detail_columns = [
        *CANONICAL_KEYS,
        "protein",
        "variant",
        *context_keys,
        "model",
        "assay_id_abundance",
        "score_abundance",
        "assay_id_function",
        "score_function",
    ]
    detail = nonidentical.loc[:, detail_columns].rename(
        columns={
            "assay_id_abundance": "assay_id_1",
            "score_abundance": "value_1",
            "assay_id_function": "assay_id_2",
            "score_function": "value_2",
        }
    )
    detail.insert(
        detail.columns.get_loc("value_2") + 1,
        "absolute_difference",
        (detail["value_1"] - detail["value_2"]).abs(),
    )
    summary = {
        "groups_checked": len(matched),
        "exact_match_count": int(exact.sum()),
        "nonidentical_count": len(nonidentical),
        "abundance_only_count": int(merged["_merge"].eq("left_only").sum()),
        "function_only_count": int(merged["_merge"].eq("right_only").sum()),
    }
    return summary, detail.reset_index(drop=True)


def orient_scores(scores: pd.DataFrame, direction: dict[str, bool]) -> pd.DataFrame:
    """Return model columns oriented so higher values mean more impairment.

    ``direction[model]`` states whether the supplied score already has that
    orientation.  The raw model columns are retained and oriented copies use
    the ``_impairment`` suffix.
    """
    result = scores.copy()
    for model, higher_is_impaired in direction.items():
        if model not in result:
            continue
        values = pd.to_numeric(result[model], errors="coerce")
        result[f"{model}_impairment"] = values if higher_is_impaired else -values
    result.attrs = scores.attrs.copy()
    result.attrs["direction"] = direction.copy()
    return result
