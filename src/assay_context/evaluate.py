"""Metrics and leakage-controlled predictors for the frozen primary cohort."""

from __future__ import annotations

from collections.abc import Sequence

import numpy as np
import pandas as pd
from sklearn.compose import ColumnTransformer
from sklearn.impute import SimpleImputer
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import average_precision_score, roc_auc_score
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import OneHotEncoder, StandardScaler

from .splits import position_folds

MODEL_COLUMNS = ("ESM1v_single_impairment", "ESM2_650M_impairment")
C_GRID = (0.01, 0.1, 1.0, 10.0)
# B1 descriptors are frozen, auditable biochemical scales.  Hydropathy is the
# Kyte-Doolittle scale (Kyte & Doolittle, J Mol Biol 1982).  Charge is a simple
# integer convention at approximately neutral pH: D/E=-1, K/R/H=+1, all
# remaining standard amino acids=0.  Differences are mutant minus wild type.
CHARGE_CONVENTION = "D/E=-1; K/R/H=+1; all other standard residues=0"
AA_CHARGE = {
    "A": 0, "C": 0, "D": -1, "E": -1, "F": 0, "G": 0, "H": 1, "I": 0,
    "K": 1, "L": 0, "M": 0, "N": 0, "P": 0, "Q": 0, "R": 1, "S": 0,
    "T": 0, "V": 0, "W": 0, "Y": 0,
}
HYDROPATHY_SCALE_NAME = "Kyte-Doolittle"
AA_HYDROPHOBICITY = {
    "A": 1.8, "C": 2.5, "D": -3.5, "E": -3.5, "F": 2.8, "G": -0.4,
    "H": -3.2, "I": 4.5, "K": -3.9, "L": 3.8, "M": 1.9, "N": -3.5,
    "P": -1.6, "Q": -3.5, "R": -4.5, "S": -0.8, "T": -0.7, "V": 4.2,
    "W": -0.9, "Y": -1.3,
}


def evaluate_binary(y: np.ndarray, score: np.ndarray) -> dict[str, float]:
    """Compute AP/AUROC and denominators for an impairment score.

    ``score`` must already follow the frozen impairment direction (higher means
    more likely impaired).  Single-class inputs are reported as invalid rather
    than receiving a misleading AUROC/AP value.
    """
    labels = np.asarray(y).astype(int, copy=False).ravel()
    values = np.asarray(score, dtype=float).ravel()
    if labels.size != values.size:
        raise ValueError("y and score must have equal length")
    if labels.size == 0:
        raise ValueError("cannot evaluate an empty cohort")
    if not np.isin(labels, (0, 1)).all():
        raise ValueError("y must contain only binary 0/1 labels")
    if not np.isfinite(values).all():
        raise ValueError("score contains non-finite values")
    positives = int(labels.sum())
    negatives = int(len(labels) - positives)
    prevalence = positives / len(labels)
    valid = positives > 0 and negatives > 0
    if valid:
        # These vectorized forms avoid repeatedly constructing sklearn metric
        # objects across 2,000 clustered draws.  Fall back to sklearn only for
        # ties, whose threshold semantics need explicit averaging.
        order = np.argsort(-values, kind="mergesort")
        sorted_scores = values[order]
        sorted_labels = labels[order]
        if np.unique(sorted_scores).size == len(sorted_scores):
            cumulative = np.cumsum(sorted_labels)
            ap = float((cumulative[sorted_labels == 1] / (np.flatnonzero(sorted_labels == 1) + 1)).sum() / positives)
            ranks = np.empty(len(labels), dtype=float)
            ranks[order] = np.arange(len(labels), 0, -1)
            auroc = float((ranks[labels == 1].sum() - positives * (positives + 1) / 2) / (positives * negatives))
        else:
            ap = float(average_precision_score(labels, values))
            auroc = float(roc_auc_score(labels, values))
        lift = float(ap / prevalence) if prevalence else np.nan
    else:
        ap = auroc = lift = np.nan
    return {
        "n_variants": float(len(labels)),
        "n_positives": float(positives),
        "n_negatives": float(negatives),
        "positives": float(positives),
        "negatives": float(negatives),
        "prevalence": float(prevalence),
        "ap": ap,
        "ap_prevalence": lift,
        "ap_lift": lift,
        "auroc": auroc,
        "valid": bool(valid),
    }


def _required_primary(frame: pd.DataFrame) -> pd.DataFrame:
    required = {"primary_eligible", "primary_endpoint", "position"}
    missing = sorted(required - set(frame.columns))
    if missing:
        raise ValueError(f"benchmark frame is missing required columns: {', '.join(missing)}")
    eligible = frame.loc[frame["primary_eligible"].astype(bool)].copy()
    if not eligible["primary_endpoint"].isin(("function-impaired", "function-preserved")).all():
        raise ValueError("primary_eligible rows must have a frozen binary primary_endpoint")
    return eligible


def assert_common_cohort(frame: pd.DataFrame, score_columns: Sequence[str] = MODEL_COLUMNS) -> None:
    """Fail if model-covered canonical variant cohorts differ."""
    id_columns = [column for column in ("protein", "sequence_hash", "wt", "position", "mutant", "variant") if column in frame]
    if not id_columns:
        id_columns = ["row_id"] if "row_id" in frame else []
    if not id_columns:
        raise ValueError("common-cohort check requires canonical variant columns")
    covered = []
    for column in score_columns:
        if column not in frame:
            raise ValueError(f"missing frozen score column: {column}")
        covered.append(set(map(tuple, frame.loc[frame[column].notna(), id_columns].astype(str).to_numpy())))
    if covered and any(values != covered[0] for values in covered[1:]):
        raise ValueError("canonical IDs differ across model score cohorts")


def evaluate_primary(
    frame: pd.DataFrame,
    score_columns: Sequence[str] = MODEL_COLUMNS,
    *,
    seed: int = 2026,
    bootstrap_draws: int = 2000,
) -> pd.DataFrame:
    """Evaluate frozen direct scores per protein, with paired position bootstrap.

    This is direct zero-shot evaluation: no split is applied and results are
    explicitly labelled as measured-position bootstrap intervals, not OOF.
    """
    eligible = _required_primary(frame)
    assert_common_cohort(eligible, score_columns)
    rows: list[dict[str, object]] = []
    from .bootstrap import cluster_bootstrap

    for protein, group in eligible.groupby("protein", sort=True):
        labels = group["primary_endpoint"].eq("function-impaired").astype(int).to_numpy()
        for model in score_columns:
            result = evaluate_binary(labels, group[model].to_numpy())
            boot = cluster_bootstrap(
                group.assign(y=labels), seed=seed, n_boot=bootstrap_draws, score_columns=[model]
            )
            valid = boot.loc[boot["valid"].astype(bool), "ap"]
            ci_low = float(valid.quantile(0.025)) if len(valid) else np.nan
            ci_high = float(valid.quantile(0.975)) if len(valid) else np.nan
            row = {"protein": protein, "model": model, **result}
            row.update({"ap_ci_low": ci_low, "ap_ci_high": ci_high,
                        "bootstrap_draws": bootstrap_draws,
                        "bootstrap_invalid_fraction": 1 - (len(valid) / bootstrap_draws),
                        "bootstrap_unstable": (1 - (len(valid) / bootstrap_draws)) > 0.10,
                        "ci_interpretation": "measured residue positions only",
                        "evaluation": "direct zero-shot; not OOF"})
            rows.append(row)
    result = pd.DataFrame(rows)
    macro: list[dict[str, object]] = []
    for model, group in result.groupby("model", sort=False):
        macro.append({"protein": "macro-average (2 proteins)", "model": model,
                      "n_proteins": int(group["protein"].nunique()),
                      "ap": float(group["ap"].mean()), "prevalence": float(group["prevalence"].mean()),
                      "ap_prevalence": float(group["ap_prevalence"].mean()),
                      "mean_per_protein_ap_lift": float(group["ap_prevalence"].mean()),
                      "auroc": float(group["auroc"].mean()), "evaluation": "macro; equal protein weight; 2 proteins only; lift is mean per-protein AP/prevalence"})
    result = pd.concat([result, pd.DataFrame(macro)], ignore_index=True)
    result.attrs["paired_models"] = list(score_columns)
    result.attrs["common_cohort"] = True
    return result


def substitution_descriptors(frame: pd.DataFrame) -> pd.DataFrame:
    """Return auditable B1 charge and Kyte-Doolittle changes."""
    missing = sorted({"wt", "mutant"} - set(frame.columns))
    if missing:
        raise ValueError(f"B1 frame is missing required columns: {', '.join(missing)}")
    out = frame[["wt", "mutant"]].astype("string").copy()
    out["charge_change"] = [AA_CHARGE.get(b, np.nan) - AA_CHARGE.get(a, np.nan) for a, b in zip(out["wt"], out["mutant"], strict=True)]
    out["hydrophobicity_change"] = [AA_HYDROPHOBICITY.get(b, np.nan) - AA_HYDROPHOBICITY.get(a, np.nan) for a, b in zip(out["wt"], out["mutant"], strict=True)]
    return out


def _feature_frame(frame: pd.DataFrame, feature_set: str) -> pd.DataFrame:
    if feature_set == "B0":
        return pd.DataFrame(index=frame.index)
    required = {"wt", "mutant"}
    if feature_set == "B2":
        required = required | set(MODEL_COLUMNS[:1])
    missing = sorted(required - set(frame.columns))
    if missing:
        raise ValueError(f"{feature_set} frame is missing required columns: {', '.join(missing)}")
    if feature_set == "B1":
        return substitution_descriptors(frame)
    if feature_set == "B2":
        return frame[[MODEL_COLUMNS[0]]].rename(columns={MODEL_COLUMNS[0]: "frozen_impairment_score"})
    raise ValueError("feature_set must be B0, B1, or B2")


def _classifier(feature_set: str, c_value: float) -> Pipeline:
    features = {"B1": ["wt", "mutant", "charge_change", "hydrophobicity_change"], "B2": ["frozen_impairment_score"]}.get(feature_set, [])
    categorical = [name for name in features if name in {"wt", "mutant"}]
    numeric = [name for name in features if name not in categorical]
    transformer = ColumnTransformer(
        [("categorical", Pipeline([("impute", SimpleImputer(strategy="most_frequent")), ("onehot", OneHotEncoder(handle_unknown="ignore", sparse_output=False))]), categorical),
         ("numeric", Pipeline([("impute", SimpleImputer(strategy="median")), ("scale", StandardScaler())]), numeric)],
        remainder="drop", verbose_feature_names_out=False,
    )
    return Pipeline([("preprocess", transformer), ("classifier", LogisticRegression(C=c_value, max_iter=2000, random_state=2026))])


def _usable_folds(frame: pd.DataFrame, n_splits: int) -> list[tuple[np.ndarray, np.ndarray]]:
    for candidate in (n_splits, 3):
        if candidate < 2 or frame["position"].nunique() < candidate:
            continue
        folds = position_folds(frame, candidate)
        if all(frame.iloc[train]["_y"].nunique() == 2 and frame.iloc[test]["_y"].nunique() == 2 for train, test in folds):
            return folds
    raise ValueError("insufficient class distribution across grouped-position folds")


def nested_grouped_cv(frame: pd.DataFrame, feature_set: str = "B2", *, outer_splits: int = 5) -> pd.DataFrame:
    """Fit B0/B1/B2 with grouped outer CV and grouped inner C selection.

    Every eligible row receives exactly one OOF prediction.  All imputation,
    one-hot encoding and scaling are fitted inside each training fold.
    """
    eligible = _required_primary(frame).reset_index(drop=True)
    eligible["_y"] = eligible["primary_endpoint"].eq("function-impaired").astype(int)
    folds = _usable_folds(eligible, outer_splits)
    outputs: list[pd.DataFrame] = []
    for fold_id, (train_idx, test_idx) in enumerate(folds):
        train = eligible.iloc[train_idx]
        test = eligible.iloc[test_idx]
        if feature_set == "B0":
            probability = float(train["_y"].mean())
            predictions = np.full(len(test), probability)
            selected_c = np.nan
        else:
            inner_folds = _usable_folds(train, 3)
            candidates: list[tuple[float, float]] = []
            for c_value in C_GRID:
                inner_scores = []
                for inner_train, inner_test in inner_folds:
                    model = _classifier(feature_set, c_value)
                    model.fit(_feature_frame(train.iloc[inner_train], feature_set), train.iloc[inner_train]["_y"])
                    inner_scores.append(average_precision_score(train.iloc[inner_test]["_y"], model.predict_proba(_feature_frame(train.iloc[inner_test], feature_set))[:, 1]))
                candidates.append((float(np.mean(inner_scores)), c_value))
            selected_c = max(candidates, key=lambda item: (item[0], -item[1]))[1]
            model = _classifier(feature_set, selected_c)
            model.fit(_feature_frame(train, feature_set), train["_y"])
            predictions = model.predict_proba(_feature_frame(test, feature_set))[:, 1]
        out = test[[column for column in ("protein", "variant", "position") if column in test]].copy()
        out["row_id"] = test.index.to_numpy()
        out["fold"] = fold_id
        out["y_true"] = test["_y"].to_numpy()
        out["prediction"] = predictions
        out["feature_set"] = feature_set
        out["selected_C"] = selected_c
        outputs.append(out)
    predictions = pd.concat(outputs, ignore_index=True).sort_values("row_id").reset_index(drop=True)
    if predictions["row_id"].nunique() != len(eligible) or set(predictions["row_id"]) != set(range(len(eligible))):
        raise AssertionError("OOF coverage is not exactly once per eligible row")
    predictions.attrs["metrics"] = evaluate_binary(predictions["y_true"].to_numpy(), predictions["prediction"].to_numpy())
    predictions.attrs["split"] = "outer grouped-position CV; nested grouped-position inner selection"
    return predictions


def leave_one_protein_out(frame: pd.DataFrame, feature_set: str = "B2") -> pd.DataFrame:
    """Evaluate only the two Gate-A PASS proteins in held-out-protein setting."""
    eligible = _required_primary(frame).reset_index(drop=True)
    if "gate_a" in eligible and not eligible["gate_a"].eq("PASS").all():
        raise ValueError("leave-one-protein-out requires Gate A PASS rows only")
    proteins = sorted(eligible["protein"].unique())
    if len(proteins) != 2:
        raise ValueError("leave-one-protein-out is limited to exactly two PASS proteins")
    outputs = []
    for held_out in proteins:
        train = eligible.loc[eligible["protein"].ne(held_out)].copy()
        test = eligible.loc[eligible["protein"].eq(held_out)].copy()
        train["_y"] = train["primary_endpoint"].eq("function-impaired").astype(int)
        test["_y"] = test["primary_endpoint"].eq("function-impaired").astype(int)
        if feature_set == "B0":
            predictions = np.full(len(test), train["_y"].mean())
            selected_c = np.nan
        else:
            inner_folds = _usable_folds(train, 3)
            scores = []
            for c_value in C_GRID:
                vals = []
                for tr, va in inner_folds:
                    model = _classifier(feature_set, c_value).fit(_feature_frame(train.iloc[tr], feature_set), train.iloc[tr]["_y"])
                    vals.append(average_precision_score(train.iloc[va]["_y"], model.predict_proba(_feature_frame(train.iloc[va], feature_set))[:, 1]))
                scores.append((np.mean(vals), c_value))
            selected_c = max(scores, key=lambda item: (item[0], -item[1]))[1]
            model = _classifier(feature_set, selected_c).fit(_feature_frame(train, feature_set), train["_y"])
            predictions = model.predict_proba(_feature_frame(test, feature_set))[:, 1]
        out = test[[column for column in ("protein", "variant", "position") if column in test]].copy()
        out["held_out_protein"] = held_out
        out["fold"] = held_out
        out["y_true"] = test["_y"].to_numpy()
        out["prediction"] = predictions
        out["feature_set"] = feature_set
        out["selected_C"] = selected_c
        outputs.append(out)
    result = pd.concat(outputs, ignore_index=True)
    result.attrs["split"] = "leave-one-protein-out; two PASS proteins only; inner grouped-position CV within training protein"
    result.attrs["metrics"] = {protein: evaluate_binary(group["y_true"].to_numpy(), group["prediction"].to_numpy()) for protein, group in result.groupby("held_out_protein")}
    return result


evaluate_supervised = nested_grouped_cv
