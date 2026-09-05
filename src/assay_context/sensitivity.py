"""Structure-added, leakage-controlled models and frozen sensitivity controls."""

from __future__ import annotations

from collections.abc import Iterable

import numpy as np
import pandas as pd
from sklearn.compose import ColumnTransformer
from sklearn.impute import SimpleImputer
from sklearn.linear_model import LogisticRegression
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import OneHotEncoder, StandardScaler

from .bootstrap import cluster_bootstrap
from .evaluate import C_GRID, MODEL_COLUMNS
from .splits import position_folds

STRUCTURE_NUMERIC = ("accessibility", "contacts", "site_distance")
STRUCTURE_CATEGORICAL = ("secondary_structure",)


def verify_external_sha256(path: str, expected_sha256: str) -> bool:
    """Return true only when an external sensitivity source matches its freeze."""
    import hashlib
    from pathlib import Path

    source = Path(path)
    if not source.is_file():
        return False
    digest = hashlib.sha256()
    with source.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest().lower() == str(expected_sha256).strip().lower()


def _ids(frame: pd.DataFrame, available: set[str] | None = None) -> list[str]:
    columns = [c for c in ("protein", "protein_id", "sequence_hash", "wt", "position", "mutant", "variant") if c in frame and (available is None or c in available)]
    if not columns:
        raise ValueError("structure cohort needs canonical variant identifiers")
    return columns


def _id_tuples(frame: pd.DataFrame, columns: list[str] | None = None) -> set[tuple[str, ...]]:
    columns = columns or _ids(frame)
    return set(map(tuple, frame[columns].astype("string").fillna("<NA>").to_numpy()))


def assert_common_variant_ids(*frames: pd.DataFrame) -> None:
    """Require exact canonical IDs across model/sensitivity cohorts."""
    if len(frames) < 2:
        raise ValueError("common-cohort check requires at least two frames")
    common_columns = set(frames[0].columns)
    for frame in frames[1:]:
        common_columns &= set(frame.columns)
    columns = _ids(frames[0], common_columns)
    if not columns:
        raise ValueError("common-cohort check requires shared canonical identifiers")
    first = _id_tuples(frames[0], columns)
    for frame in frames[1:]:
        if _id_tuples(frame, columns) != first:
            raise ValueError("canonical variant cohort differs across models")


def build_structure_features(frame: pd.DataFrame) -> pd.DataFrame:
    """Build B3-compatible features and explicit missingness indicators."""
    missing = sorted(set(STRUCTURE_NUMERIC + STRUCTURE_CATEGORICAL) - set(frame.columns))
    if missing:
        raise ValueError(f"structure frame is missing features: {', '.join(missing)}")
    output = frame.loc[:, list(STRUCTURE_NUMERIC) + list(STRUCTURE_CATEGORICAL)].copy()
    for column in STRUCTURE_NUMERIC:
        output[column] = pd.to_numeric(output[column], errors="coerce")
        output[f"{column}_missing"] = output[column].isna()
    for column in STRUCTURE_CATEGORICAL:
        output[column] = output[column].astype("string").fillna("missing")
        output[f"{column}_missing"] = output[column].eq("missing")
    return output


def structure_complete_case(frame: pd.DataFrame) -> pd.DataFrame:
    """Return rows with all measured structural predictors available."""
    required = set(STRUCTURE_NUMERIC + STRUCTURE_CATEGORICAL)
    missing = sorted(required - set(frame.columns))
    if missing:
        raise ValueError(f"structure frame is missing features: {', '.join(missing)}")
    mask = frame[list(STRUCTURE_NUMERIC)].apply(pd.to_numeric, errors="coerce").notna().all(axis=1)
    mask &= frame[list(STRUCTURE_CATEGORICAL)].notna().all(axis=1)
    return frame.loc[mask].copy()


def _feature_frame(frame: pd.DataFrame, feature_set: str) -> pd.DataFrame:
    structure = build_structure_features(frame)
    if feature_set == "B3":
        return structure
    if feature_set == "B4":
        if MODEL_COLUMNS[0] not in frame:
            raise ValueError("B4 frame is missing frozen ESM1v impairment score")
        structure["frozen_impairment_score"] = pd.to_numeric(frame[MODEL_COLUMNS[0]], errors="coerce")
        structure["frozen_impairment_score_missing"] = structure["frozen_impairment_score"].isna()
        return structure
    raise ValueError("feature_set must be B3 or B4")


def _model(feature_set: str, c_value: float) -> Pipeline:
    numeric = [*STRUCTURE_NUMERIC, *(f"{c}_missing" for c in STRUCTURE_NUMERIC), *(f"{c}_missing" for c in STRUCTURE_CATEGORICAL)]
    if feature_set == "B4":
        numeric += ["frozen_impairment_score", "frozen_impairment_score_missing"]
    categorical = list(STRUCTURE_CATEGORICAL)
    transformer = ColumnTransformer(
        [
            ("numeric", Pipeline([("impute", SimpleImputer(strategy="median")), ("scale", StandardScaler())]), numeric),
            ("categorical", Pipeline([("impute", SimpleImputer(strategy="most_frequent")), ("onehot", OneHotEncoder(handle_unknown="ignore", sparse_output=False))]), categorical),
        ],
        remainder="drop", verbose_feature_names_out=False,
    )
    return Pipeline([("preprocess", transformer), ("classifier", LogisticRegression(C=c_value, max_iter=2000, random_state=2026))])


def _validate_folds(frame: pd.DataFrame, fold_column: str) -> list[tuple[np.ndarray, np.ndarray]]:
    folds = []
    values = pd.to_numeric(frame[fold_column], errors="raise").astype(int)
    if values.isna().any() or values.nunique() < 2:
        raise ValueError("frozen fold assignment must contain at least two folds")
    for fold in sorted(values.unique()):
        test = np.flatnonzero(values.to_numpy() == fold)
        train = np.flatnonzero(values.to_numpy() != fold)
        group_columns = ["position"]
        if "protein" in frame:
            group_columns.insert(0, "protein")
        train_groups = set(map(tuple, frame.iloc[train][group_columns].astype("string").to_numpy()))
        test_groups = set(map(tuple, frame.iloc[test][group_columns].astype("string").to_numpy()))
        if train_groups & test_groups:
            raise ValueError("frozen folds leak residue positions")
        folds.append((train, test))
    if len(np.concatenate([test for _, test in folds])) != len(frame) or len(set(np.concatenate([test for _, test in folds]))) != len(frame):
        raise ValueError("frozen folds do not cover every row exactly once")
    return folds


def compare_structure_models(frame: pd.DataFrame, *, fold_column: str = "fold", feature_sets: Iterable[str] = ("B3", "B4")) -> pd.DataFrame:
    """Generate structure-only and score-plus-structure OOF predictions.

    ``fold_column`` is the exact outer fold assignment generated by Task 4;
    fitting (including imputation, scaling and one-hot encoding) happens only
    after the test rows have been separated.
    """
    required = {"position", "primary_endpoint", fold_column}
    missing = sorted(required - set(frame.columns))
    if missing:
        raise ValueError(f"structure benchmark frame is missing: {', '.join(missing)}")
    eligibility = frame["primary_eligible"].astype(bool) if "primary_eligible" in frame else pd.Series(True, index=frame.index)
    eligible = frame.loc[eligibility].reset_index(drop=True)
    eligible = eligible.loc[eligible["primary_endpoint"].isin(("function-impaired", "function-preserved"))].reset_index(drop=True)
    folds = _validate_folds(eligible, fold_column)
    outputs: list[pd.DataFrame] = []
    for feature_set in feature_sets:
        if feature_set not in {"B3", "B4"}:
            raise ValueError("feature_sets may contain only B3 and B4")
        for fold_id, (train_idx, test_idx) in enumerate(folds):
            train, test = eligible.iloc[train_idx], eligible.iloc[test_idx]
            y_train = train["primary_endpoint"].eq("function-impaired").astype(int)
            if y_train.nunique() < 2:
                # A fold with one training class is an explicit missing result,
                # not a leakage-prone fallback model.
                predictions = np.full(len(test), np.nan)
                selected_c = np.nan
            else:
                candidates: list[tuple[float, float]] = []
                inner_folds = position_folds(train, min(3, train["position"].nunique())) if train["position"].nunique() >= 2 else []
                for c_value in C_GRID:
                    inner_scores = []
                    for inner_train, inner_test in inner_folds:
                        y_inner = y_train.iloc[inner_train]
                        if y_inner.nunique() < 2:
                            continue
                        inner_model = _model(feature_set, c_value)
                        inner_model.fit(_feature_frame(train.iloc[inner_train], feature_set), y_inner)
                        probabilities = inner_model.predict_proba(_feature_frame(train.iloc[inner_test], feature_set))[:, 1]
                        y_test = y_train.iloc[inner_test]
                        if y_test.nunique() == 2:
                            from sklearn.metrics import average_precision_score

                            inner_scores.append(float(average_precision_score(y_test, probabilities)))
                    # If a tiny training set cannot provide a two-class inner
                    # test fold, select C by training-only AP as an explicit,
                    # deterministic fallback (never by outer-test labels).
                    if inner_scores:
                        candidates.append((float(np.mean(inner_scores)), c_value))
                    else:
                        model = _model(feature_set, c_value)
                        model.fit(_feature_frame(train, feature_set), y_train)
                        candidates.append((float(model.score(_feature_frame(train, feature_set), y_train)), c_value))
                selected_c = max(candidates, key=lambda item: (item[0], -item[1]))[1]
                model = _model(feature_set, selected_c)
                model.fit(_feature_frame(train, feature_set), y_train)
                predictions = model.predict_proba(_feature_frame(test, feature_set))[:, 1]
            columns = [c for c in ("protein", "protein_id", "sequence_hash", "variant", "wt", "position", "mutant") if c in test]
            output = test[columns].copy()
            output["row_id"] = test.index.to_numpy()
            output["fold"] = pd.to_numeric(test[fold_column], errors="raise").astype(int).to_numpy()
            output["y_true"] = test["primary_endpoint"].eq("function-impaired").astype(int).to_numpy()
            output["prediction"] = predictions
            output["feature_set"] = feature_set
            output["selected_C"] = selected_c
            outputs.append(output)
    result = pd.concat(outputs, ignore_index=True)
    result.attrs["split"] = "Task 4 outer grouped-position folds reused exactly; preprocessing fit within each training fold"
    result.attrs["feature_sets"] = tuple(feature_sets)
    return result


def assign_proximity(frame: pd.DataFrame, threshold: float) -> pd.Series:
    """Return a descriptive proximity flag from minimum heavy-atom distance."""
    if "site_distance" not in frame:
        raise ValueError("site_distance is required for proximity")
    return pd.to_numeric(frame["site_distance"], errors="coerce").le(float(threshold)) & frame["site_distance"].notna()


def sensitivity_statuses(*, function_source_available: bool = True, structure_available: bool = True) -> pd.DataFrame:
    """Enumerate every frozen control, including non-applicable controls."""
    rows = [
        ("function_score_source", "PASS" if function_source_available else "FAILED", "external function-assay score source hash validated before use" if function_source_available else "external source missing or hash mismatch"),
        ("label_confidence", "NOT_APPLICABLE", "quantitative confidence-regime comparison requires a frame; do not infer PASS from retained columns"),
        ("coverage", "NOT_APPLICABLE", "quantitative common-versus-available counts require a frame; use task5_sensitivity_results.csv"),
        ("structure_missingness", "NOT_APPLICABLE", "quantitative complete-case versus imputed comparison requires a frame; missing residues are retained as NaN" if structure_available else "no validated structure features available"),
        ("study_system", "NOT_APPLICABLE", "same-system/study comparison is confounded unless an unconfounded frame is supplied"),
        ("proximity_4A", "PASS" if structure_available else "NOT_APPLICABLE", "descriptive minimum heavy-atom distance threshold requires denominators"),
        ("proximity_5A", "PASS" if structure_available else "NOT_APPLICABLE", "predeclared descriptive threshold requires denominators"),
        ("proximity_8A", "PASS" if structure_available else "NOT_APPLICABLE", "sensitivity threshold requires denominators; not allostery"),
        ("threshold_variations", "NOT_APPLICABLE", "no post-hoc alternatives were frozen before Task 4"),
        ("permuted_labels", "NOT_APPLICABLE", "quantitative Task 4 permutation summary and verified hash are required"),
        ("protein_opposite_or_no_improvement", "NOT_APPLICABLE", "paired effect and confidence interval are required"),
    ]
    return pd.DataFrame(rows, columns=["control", "status", "reason"])


def paired_prediction_bootstrap(frame: pd.DataFrame, *, setting: str, seed: int = 2026, n_boot: int = 2000) -> pd.DataFrame:
    """Bootstrap B2/B3/B4 predictions with one shared position draw.

    The input is long by feature set; canonical IDs and labels must be exactly
    repeated for each model.  One ``cluster_bootstrap`` call per protein keeps
    the sampled position sequence paired across all three predictors.
    """
    required = {"protein", "variant", "position", "y_true", "feature_set", "prediction"}
    missing = sorted(required - set(frame.columns))
    if missing:
        raise ValueError(f"prediction bootstrap frame is missing: {', '.join(missing)}")
    models = ("B2", "B3", "B4")
    model_ids = {
        model: _id_tuples(frame.loc[frame["feature_set"].eq(model)], ["protein", "variant", "position"])
        for model in models
    }
    if any(model_ids[model] != model_ids["B2"] for model in models[1:]):
        raise ValueError("prediction bootstrap requires exact common canonical IDs across B2/B3/B4")
    wide = frame.loc[frame["feature_set"].isin(models)].pivot_table(
        index=["protein", "variant", "position", "y_true"], columns="feature_set", values="prediction", aggfunc="first"
    ).reset_index()
    if set(wide.columns) & set(models) != set(models):
        raise ValueError("prediction bootstrap requires B2, B3, and B4")
    model_rows = len(frame.loc[frame["feature_set"].isin(models)])
    if model_rows != 3 * len(wide):
        raise ValueError("prediction bootstrap has duplicate or incomplete canonical model rows")
    outputs = []
    for protein, group in wide.groupby("protein", sort=True):
        draws = cluster_bootstrap(group.assign(y=group["y_true"]), seed, n_boot, score_columns=list(models))
        draws.insert(0, "protein", protein)
        draws.insert(1, "setting", setting)
        outputs.append(draws)
    return pd.concat(outputs, ignore_index=True)


def validate_sensitivity_evidence(frame: pd.DataFrame) -> None:
    """Fail closed if a PASS row lacks quantitative evidence."""
    required = {"control", "status", "reason", "n_variants", "n_positions", "metric", "effect", "evidence"}
    missing = sorted(required - set(frame.columns))
    if missing:
        raise ValueError(f"sensitivity evidence is missing numeric evidence fields: {', '.join(missing)}")
    if not frame["status"].isin(("PASS", "NOT_APPLICABLE", "FAILED")).all():
        raise ValueError("unknown sensitivity status")
    passed = frame.loc[frame["status"].eq("PASS")]
    if passed.empty:
        return
    for column in ("n_variants", "n_positions", "metric", "effect"):
        if passed[column].isna().any():
            raise ValueError("PASS sensitivity row lacks numeric evidence")
    if passed["evidence"].astype("string").str.strip().eq("").any():
        raise ValueError("PASS sensitivity row lacks numeric evidence")


def assert_complete_b4_b2_contrasts(frame: pd.DataFrame) -> None:
    """Require one paired B4-minus-B2 contrast for every plot stratum."""
    required = {("within-protein", "CYP2C9"), ("within-protein", "PTEN"),
                ("held-out protein", "CYP2C9"), ("held-out protein", "PTEN")}
    columns = {"setting", "protein", "model"}
    missing = sorted(columns - set(frame.columns))
    if missing:
        raise ValueError(f"paired contrast metadata is missing: {', '.join(missing)}")
    actual = set(map(tuple, frame.loc[frame["model"].eq("B4-B2"), ["setting", "protein"]].drop_duplicates().to_numpy()))
    if actual != required:
        raise ValueError(f"paired B4-B2 contrast strata differ from required four: {sorted(actual)}")
