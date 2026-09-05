"""Author-calibrated phenotype adapters used for the frozen paired cohort."""

from __future__ import annotations

import math
from collections.abc import Iterable

import pandas as pd

# Exact 0.975 two-sided t critical values for the finite df values present in
# the GCK source tables. Unexpected values fail explicitly at the freeze gate.
_T975 = {
    1: 12.706, 2: 4.303, 3: 3.182, 4: 2.776, 5: 2.571, 6: 2.447,
    7: 2.365, 8: 2.306, 9: 2.262, 10: 2.228, 11: 2.201, 12: 2.179,
    13: 2.160, 14: 2.145, 15: 2.131, 16: 2.120, 17: 2.110, 18: 2.101,
    19: 2.093, 20: 2.086,
}


def label_interval(low: float, high: float, threshold: float) -> str:
    """Label an interval as low/preserved only when wholly separated.

    Equality at a threshold is intentionally uncertain.  This makes boundary
    behavior explicit and avoids turning a rounded confidence limit into a
    falsely certain class.
    """
    values = (low, high, threshold)
    try:
        converted = tuple(float(value) for value in values)
    except (TypeError, ValueError):
        return "uncertain"
    if any(not math.isfinite(value) for value in converted):
        return "uncertain"
    low_value, high_value, threshold_value = converted
    if low_value > high_value:
        raise ValueError("Interval lower bound must not exceed upper bound")
    if high_value < threshold_value:
        return "low"
    if low_value > threshold_value:
        return "preserved"
    return "uncertain"


def gate_a_status(frame: pd.DataFrame) -> str:
    """Apply all frozen Gate A minimums to a labeled primary-endpoint frame."""
    required = {"position", "primary_endpoint"}
    missing = sorted(required - set(frame.columns))
    if missing:
        raise ValueError(f"Gate A frame missing required columns: {', '.join(missing)}")
    if len(frame) < 300 or frame["position"].nunique() < 50:
        return "REDUCE_SCOPE"
    class_rows = {
        label: frame.loc[frame["primary_endpoint"].eq(label)]
        for label in ("function-impaired", "function-preserved")
    }
    if any(len(rows) < 50 or rows["position"].nunique() < 10 for rows in class_rows.values()):
        return "REDUCE_SCOPE"
    return "PASS"


def primary_binary_eligibility(frame: pd.DataFrame, gate_a: str) -> pd.Series:
    """Return rows eligible for the frozen benchmark only after Gate A passes."""
    if "primary_endpoint" not in frame:
        raise ValueError("Benchmark frame missing primary_endpoint column")
    return frame["primary_endpoint"].isin(("function-impaired", "function-preserved")) & (gate_a == "PASS")


def _require_columns(frame: pd.DataFrame, required: Iterable[str], adapter: str) -> None:
    missing = sorted(set(required) - set(frame.columns))
    if missing:
        raise ValueError(f"{adapter} source missing required columns: {', '.join(missing)}")


def _variant_column(frame: pd.DataFrame, adapter: str) -> str:
    for column in ("variant", "Variant (one letter)"):
        if column in frame.columns:
            return column
    raise ValueError(f"{adapter} source missing variant column")


def _prepare(frame: pd.DataFrame, adapter: str, source_class: str | None = None) -> pd.DataFrame:
    column = _variant_column(frame, adapter)
    if source_class is not None:
        _require_columns(frame, (source_class,), adapter)
        frame = frame.loc[frame[source_class].astype("string").str.lower().eq("missense")].copy()
    else:
        frame = frame.copy()
    frame = frame.rename(columns={column: "variant"})
    frame["variant"] = frame["variant"].astype("string").str.strip().str.upper()
    if frame["variant"].isna().any() or frame["variant"].eq("").any():
        raise ValueError(f"{adapter} source contains missing variant values")
    if frame["variant"].duplicated().any():
        duplicates = frame.loc[frame["variant"].duplicated(keep=False), "variant"].drop_duplicates().tolist()
        raise ValueError(f"{adapter} source contains duplicate variants: {duplicates}")
    return frame.reset_index(drop=True)


def _class_label(value: object, preserved: str, low: str, *, high_is_other: bool = True) -> str:
    if value is None or pd.isna(value):
        return "uncertain"
    label = str(value).strip().lower()
    if label == preserved:
        return "preserved"
    if label == low:
        return "low" if low == "low" else "impaired"
    if label.startswith("possibly_") or label in {"unknown", "missing", "na"}:
        return "uncertain"
    if high_is_other and label in {"increased", "high", "possibly_high"}:
        return "other"
    return "uncertain"


def _function_class_label(value: object) -> str:
    if value is None or pd.isna(value):
        return "uncertain"
    label = str(value).strip().lower()
    if label == "wt-like":
        return "preserved"
    if label in {"decreased", "nonsense-like"}:
        return "impaired"
    if label.startswith("possibly_"):
        return "uncertain"
    if label in {"increased", "high"}:
        return "other"
    return "uncertain"


def _abundance_class_label(value: object) -> str:
    if value is None or pd.isna(value):
        return "uncertain"
    label = str(value).strip().lower()
    if label == "wt-like":
        return "preserved"
    if label in {"decreased", "nonsense-like"}:
        return "low"
    if label.startswith("possibly_"):
        return "uncertain"
    if label in {"increased", "high"}:
        return "other"
    return "uncertain"


def _class_result(frame: pd.DataFrame, abundance_column: str, function_column: str, adapter: str) -> pd.DataFrame:
    result = pd.DataFrame({"variant": frame["variant"]})
    result["abundance_label"] = frame[abundance_column].map(_abundance_class_label)
    result["function_label"] = frame[function_column].map(_function_class_label)
    result["uncertainty_kind"] = "point-estimate author category"
    result["source_adapter"] = adapter
    return result


def adapt_cyp2c9(source: pd.DataFrame) -> pd.DataFrame:
    """Adapt Amorosi CYP2C9 classes; increased/possibly classes stay out."""
    frame = _prepare(source, "CYP2C9", "class" if "class" in source.columns else None)
    _require_columns(frame, ("abundance_class", "activity_class"), "CYP2C9")
    return _class_result(frame, "abundance_class", "activity_class", "cyp2c9_author_classes")


def adapt_pten_abundance(source: pd.DataFrame) -> pd.DataFrame:
    frame = _prepare(source, "PTEN abundance", "class" if "class" in source.columns else None)
    _require_columns(frame, ("abundance_class",), "PTEN abundance")
    result = pd.DataFrame({"variant": frame["variant"]})
    result["abundance_label"] = frame["abundance_class"].map(
        lambda value: _class_label(value, "wt-like", "low")
    )
    result["uncertainty_kind"] = "point-estimate author category"
    result["source_adapter"] = "pten_composite_abundance_author_classes"
    return result


def _as_bool(value: object) -> bool:
    if isinstance(value, bool):
        return value
    if value is None or pd.isna(value):
        return False
    return str(value).strip().lower() in {"true", "t", "1", "yes"}


def adapt_pten_function(source: pd.DataFrame) -> pd.DataFrame:
    frame = _prepare(source, "PTEN function", "Type" if "Type" in source.columns else None)
    _require_columns(frame, ("Cum_score", "High_conf"), "PTEN function")
    score = pd.to_numeric(frame["Cum_score"], errors="coerce")
    result = pd.DataFrame({"variant": frame["variant"]})
    labels: list[str] = []
    for value, high_conf in zip(score, frame["High_conf"], strict=True):
        if not _as_bool(high_conf) or pd.isna(value):
            labels.append("uncertain")
        elif value <= -1.11:
            labels.append("impaired")
        elif value <= 0.89:
            labels.append("preserved")
        else:
            labels.append("other")
    result["function_label"] = labels
    result["function_score"] = score
    result["uncertainty_kind"] = "point-estimate high-confidence author score"
    result["source_adapter"] = "pten_mighell_high_conf_cutoffs"
    return result


def _t_critical(df: object) -> float | None:
    if df is None:
        return None
    try:
        if bool(pd.isna(df)):
            return None
    except (TypeError, ValueError):
        raise ValueError(f"unsupported GCK degrees of freedom: {df!r}") from None
    try:
        value = float(df)
    except (TypeError, ValueError):
        raise ValueError(f"unsupported GCK degrees of freedom: {df!r}") from None
    if value <= 0 or not value.is_integer():
        raise ValueError(f"unsupported GCK degrees of freedom: {df!r}")
    integer = int(value)
    if integer in _T975:
        return _T975[integer]
    raise ValueError(f"unsupported GCK degrees of freedom: {df!r}")


def _ci(score: object, se: object, df: object) -> tuple[float | None, float | None]:
    score_value = pd.to_numeric(pd.Series([score]), errors="coerce").iloc[0]
    se_value = pd.to_numeric(pd.Series([se]), errors="coerce").iloc[0]
    critical = _t_critical(df)
    if pd.isna(score_value) or pd.isna(se_value) or critical is None:
        return None, None
    margin = critical * abs(float(se_value))
    return float(score_value) - margin, float(score_value) + margin


def _gck_label(low: float | None, high: float | None, role: str) -> str:
    if low is None or high is None:
        return "uncertain"
    if role == "abundance":
        return label_interval(low, high, 0.6)
    # Equality to either source threshold is uncertain, not a certain class.
    if high < 0.66:
        return "impaired"
    if low > 0.66 and high < 1.18:
        return "preserved"
    if low > 1.18:
        return "other"
    return "uncertain"


def adapt_gck(source: pd.DataFrame, *, role: str) -> pd.DataFrame:
    if role not in {"abundance", "function"}:
        raise ValueError("GCK role must be abundance or function")
    frame = _prepare(source, f"GCK {role}")
    score_column = "abundance_score" if role == "abundance" else "activity_score"
    se_column = f"{score_column}_se"
    _require_columns(frame, (score_column, se_column, "df"), f"GCK {role}")
    intervals = [_ci(score, se, df) for score, se, df in zip(frame[score_column], frame[se_column], frame["df"], strict=True)]
    result = pd.DataFrame({"variant": frame["variant"]})
    result[f"{role}_score"] = pd.to_numeric(frame[score_column], errors="coerce")
    result[f"{role}_lower_ci"] = [low for low, _ in intervals]
    result[f"{role}_upper_ci"] = [high for _, high in intervals]
    result[f"{role}_label"] = [_gck_label(low, high, role) for low, high in intervals]
    result["uncertainty_kind"] = "95% CI t(df), equality uncertain"
    result["source_adapter"] = f"gck_{role}_ci_thresholds"
    return result


def adapt_vkor(source: pd.DataFrame) -> pd.DataFrame:
    frame = _prepare(source, "VKORC1", "class" if "class" in source.columns else None)
    _require_columns(frame, ("abundance_class", "activity_class"), "VKORC1")
    result = _class_result(frame, "abundance_class", "activity_class", "vkor_author_classes")
    result["abundance_label"] = frame["abundance_class"].map(
        lambda value: _class_label(value, "wt-like", "low")
    )
    result["function_label"] = frame["activity_class"].map(
        lambda value: _class_label(value, "wt-like", "low")
    ).replace({"low": "impaired"})
    return result
