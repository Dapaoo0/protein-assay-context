"""Descriptive phenotype atlas tables and position-cluster intervals."""

from __future__ import annotations

import numpy as np
import pandas as pd

PHENOTYPES = (
    "abundance-preserved_function-impaired",
    "abundance-preserved_function-preserved",
    "abundance-low_function-impaired",
    "abundance-low_function-preserved",
    "uncertain-or-other",
)


def _category(abundance: object, function: object) -> str:
    a, f = str(abundance).lower(), str(function).lower()
    if a == "preserved" and f == "impaired":
        return PHENOTYPES[0]
    if a == "preserved" and f == "preserved":
        return PHENOTYPES[1]
    if a == "low" and f == "impaired":
        return PHENOTYPES[2]
    if a == "low" and f == "preserved":
        return PHENOTYPES[3]
    return PHENOTYPES[4]


def _position_bootstrap(
    frame: pd.DataFrame, categories: pd.Series, *, seed: int, draws: int
) -> dict[str, tuple[float, float]]:
    """Bootstrap fractions by sampling residue positions as intact clusters."""
    rng = np.random.default_rng(seed)
    positions = frame["position"].drop_duplicates().to_numpy()
    groups = [frame.index[frame["position"].eq(position)].to_numpy() for position in positions]
    fractions = {category: np.empty(draws, dtype=float) for category in PHENOTYPES}
    for draw in range(draws):
        selected = rng.integers(0, len(groups), size=len(groups))
        indices = np.concatenate([groups[index] for index in selected])
        values = categories.loc[indices]
        denominator = len(indices)
        for category in PHENOTYPES:
            fractions[category][draw] = float((values == category).sum()) / denominator
    return {
        category: (float(np.quantile(values, 0.025)), float(np.quantile(values, 0.975)))
        for category, values in fractions.items()
    }


def summarize_phenotypes(
    pairs: pd.DataFrame, *, seed: int = 2026, bootstrap_draws: int = 2000
) -> pd.DataFrame:
    """Return one row per protein and phenotype with explicit denominators."""
    required = {"protein", "position", "abundance_label", "function_label"}
    missing = sorted(required - set(pairs.columns))
    if missing:
        raise ValueError(f"pairs is missing required columns: {', '.join(missing)}")
    frame = pairs.copy().reset_index(drop=True)
    categories = frame.apply(
        lambda row: _category(row["abundance_label"], row["function_label"]), axis=1
    )
    frame["phenotype"] = categories
    rows: list[dict[str, object]] = []
    for protein, group in frame.groupby("protein", sort=True):
        group_categories = group["phenotype"]
        denominator = len(group)
        cis = _position_bootstrap(group, group_categories, seed=seed, draws=bootstrap_draws)
        for category in PHENOTYPES:
            count = int(group_categories.eq(category).sum())
            low, high = cis[category]
            rows.append(
                {
                    "protein": protein,
                    "phenotype": category,
                    "count": count,
                    "n": count,
                    "denominator": denominator,
                    "fraction": count / denominator if denominator else np.nan,
                    "position_count": int(group["position"].nunique()),
                    "fraction_ci_low": low,
                    "fraction_ci_high": high,
                    "ci_method": "95% position-bootstrap CI",
                    "bootstrap_seed": seed,
                    "bootstrap_draws": bootstrap_draws,
                }
            )
    result = pd.DataFrame(rows)
    result.attrs["denominator_definition"] = "all paired missense rows for the protein"
    return result
