"""Leakage-controlled grouped cross-validation helpers."""

from __future__ import annotations

import numpy as np
import pandas as pd
from sklearn.model_selection import GroupKFold


def position_folds(frame: pd.DataFrame, n_splits: int) -> list[tuple[np.ndarray, np.ndarray]]:
    """Return deterministic folds whose residue positions never cross a split.

    The returned arrays contain original row positions (not dataframe labels),
    which makes them safe to pass to ``iloc``.  Position is the grouping unit
    mandated by the analysis freeze; if multiple proteins are supplied, the
    protein is included in the group key to avoid cross-protein leakage at the
    same residue number.
    """
    if not isinstance(frame, pd.DataFrame):
        raise TypeError("frame must be a pandas DataFrame")
    if "position" not in frame:
        raise ValueError("frame is missing required column: position")
    if not isinstance(n_splits, (int, np.integer)) or n_splits < 2:
        raise ValueError("n_splits must be an integer >= 2")
    if len(frame) == 0:
        raise ValueError("cannot split an empty frame")
    if frame["position"].isna().any():
        raise ValueError("position contains missing values")
    if "protein" in frame and frame["protein"].nunique(dropna=False) > 1:
        groups = frame["protein"].astype("string") + "::" + frame["position"].astype("string")
    else:
        groups = frame["position"]
    n_groups = groups.nunique(dropna=False)
    if n_splits > n_groups:
        raise ValueError(f"n_splits={n_splits} exceeds {n_groups} residue-position groups")
    splitter = GroupKFold(n_splits=n_splits)
    indices = np.arange(len(frame), dtype=int)
    return [(train.astype(int), test.astype(int)) for train, test in splitter.split(indices, groups=groups)]
