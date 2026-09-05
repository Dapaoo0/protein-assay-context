"""Run the frozen direct-score and leakage-controlled predictor benchmark."""

from __future__ import annotations

import argparse
import ast
import hashlib
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from assay_context.artifacts import external_artifact_path
from assay_context.bootstrap import cluster_bootstrap, position_block_permutation
from assay_context.evaluate import (
    MODEL_COLUMNS,
    evaluate_binary,
    evaluate_primary,
    leave_one_protein_out,
    nested_grouped_cv,
)

PASS_PROTEINS = ("CYP2C9", "PTEN")
FROZEN_JOIN_KEYS = ("protein_id", "sequence_hash", "wt", "position", "mutant")
FROZEN_PROTEIN_IDS = {"CYP2C9": "CP2C9_HUMAN", "GCK": "HXK4_HUMAN", "PTEN": "PTEN_HUMAN", "VKORC1": "VKOR1_HUMAN"}


def _validated_canonical(frame: pd.DataFrame, label: str) -> pd.DataFrame:
    """Normalize and validate the frozen key before any score join."""
    required = set(FROZEN_JOIN_KEYS) | {"variant"}
    missing = sorted(required - set(frame.columns))
    if missing:
        raise ValueError(f"{label} is missing frozen canonical columns: {', '.join(missing)}")
    result = frame.copy()
    for column in ("protein_id", "sequence_hash", "wt", "mutant", "variant"):
        result[column] = result[column].astype("string").str.strip()
    result["sequence_hash"] = result["sequence_hash"].str.lower()
    result["wt"] = result["wt"].str.upper()
    result["mutant"] = result["mutant"].str.upper()
    result["variant"] = result["variant"].str.upper()
    result["position"] = pd.to_numeric(result["position"], errors="raise").astype("int64")
    if result[list(required)].isna().any().any() or result["variant"].eq("").any():
        raise ValueError(f"{label} contains missing frozen canonical values")
    expected = result["wt"] + result["position"].astype("string") + result["mutant"]
    if not result["variant"].eq(expected).all():
        raise ValueError(f"{label} has inconsistent variant display token")
    if result.duplicated(list(FROZEN_JOIN_KEYS)).any():
        raise ValueError(f"{label} contains duplicate frozen canonical keys")
    return result


def _read_primary(cohort_path: Path, atlas_path: Path) -> pd.DataFrame:
    cohort_raw = pd.read_csv(cohort_path)
    if "protein_id" not in cohort_raw:
        if "protein" not in cohort_raw:
            raise ValueError("cohort is missing protein_id and protein mapping column")
        cohort_raw["protein_id"] = cohort_raw["protein"].map(FROZEN_PROTEIN_IDS)
    cohort = _validated_canonical(cohort_raw, "cohort")
    atlas = _validated_canonical(pd.read_csv(atlas_path), "atlas")
    score_columns = [*FROZEN_JOIN_KEYS, "variant", "protein", *MODEL_COLUMNS]
    missing = sorted(set(score_columns) - set(atlas.columns))
    if missing:
        raise ValueError(f"atlas is missing score columns: {', '.join(missing)}")
    joined = cohort.merge(atlas[score_columns], on=list(FROZEN_JOIN_KEYS), how="left", suffixes=("_cohort", "_atlas"), validate="one_to_one")
    if "variant_atlas" in joined and not joined["variant_atlas"].eq(joined["variant_cohort"]).all():
        raise ValueError("atlas variant display token mismatch for frozen canonical key")
    if "protein_atlas" in joined and not joined["protein_atlas"].eq(joined["protein_cohort"]).all():
        raise ValueError("atlas protein display mismatch for frozen canonical key")
    joined = joined.rename(columns={"protein_cohort": "protein", "variant_cohort": "variant"})
    # A same-token row that misses the full key is an attempted key downgrade
    # (for example, a different sequence hash), not ordinary score missingness.
    atlas_tokens = set(zip(atlas["protein"].astype(str), atlas["variant"].astype(str), strict=True))
    unmatched_tokens = set(zip(cohort["protein"].astype(str), cohort["variant"].astype(str), strict=True)) & atlas_tokens
    matched_tokens = set(zip(joined.loc[joined[MODEL_COLUMNS[0]].notna(), "protein"].astype(str), joined.loc[joined[MODEL_COLUMNS[0]].notna(), "variant"].astype(str), strict=True))
    if unmatched_tokens - matched_tokens:
        raise ValueError("canonical key mismatch between cohort and atlas (same variant token, different frozen key)")
    # The direct primary cohort is the frozen two-protein subset, not all rows
    # carrying author labels or all descriptive proteins.
    eligible = joined.loc[joined["primary_eligible"].astype(bool) & joined["protein"].isin(PASS_PROTEINS)].copy()
    for column in MODEL_COLUMNS:
        if eligible[column].isna().any():
            raise ValueError(f"primary score coverage is incomplete for {column}")
    return eligible


def _direct_bootstrap(eligible: pd.DataFrame, seed: int, draws: int) -> pd.DataFrame:
    chunks = []
    for protein, group in eligible.groupby("protein", sort=True):
        labels = group["primary_endpoint"].eq("function-impaired").astype(int).to_numpy()
        boot = cluster_bootstrap(group.assign(y=labels), seed, draws, score_columns=MODEL_COLUMNS)
        boot.insert(0, "protein", protein)
        chunks.append(boot)
    result = pd.concat(chunks, ignore_index=True)
    paired = result.pivot(index=["protein", "draw"], columns="model", values="ap").reset_index()
    paired.columns.name = None
    paired["paired_ap_contrast_esm2_minus_esm1"] = paired[MODEL_COLUMNS[1]] - paired[MODEL_COLUMNS[0]]
    return result.merge(paired[["protein", "draw", "paired_ap_contrast_esm2_minus_esm1"]], on=["protein", "draw"], how="left")


def _bootstrap_summary(direct: pd.DataFrame, bootstrap: pd.DataFrame, output_path: Path, seed: int, draws: int) -> pd.DataFrame:
    """Write draw-level rows externally and retain only an auditable summary."""
    output_path.parent.mkdir(parents=True, exist_ok=True)
    bootstrap.to_csv(output_path, index=False)
    digest = hashlib.sha256(output_path.read_bytes()).hexdigest()
    rows: list[dict[str, object]] = []
    for row in direct.loc[direct["protein"].isin(PASS_PROTEINS)].itertuples(index=False):
        sample = bootstrap.loc[(bootstrap["protein"].eq(row.protein)) & (bootstrap["model"].eq(row.model))]
        valid = sample.loc[sample["valid"].astype(bool), "ap"]
        rows.append({"protein": row.protein, "model": row.model, "n_variants": int(row.n_variants), "n_positions": int(sample["sampled_positions"].iloc[0].__len__()), "prevalence": row.prevalence, "ap": row.ap, "ap_ci_low": valid.quantile(0.025), "ap_ci_high": valid.quantile(0.975), "invalid_fraction": 1 - len(valid) / draws, "unstable": (1 - len(valid) / draws) > 0.10, "seed": seed, "n_boot": draws, "ci_interpretation": "measured residue positions only", "full_draw_output": str(output_path), "full_draw_sha256": digest, "full_draw_rows": len(bootstrap)})
    paired = bootstrap.groupby(["protein", "draw"], as_index=False)["paired_ap_contrast_esm2_minus_esm1"].first()
    for protein, group in paired.groupby("protein", sort=True):
        all_values = group["paired_ap_contrast_esm2_minus_esm1"]
        values = all_values.dropna()
        invalid_fraction = float(all_values.isna().sum() / len(all_values))
        rows.append({"protein": protein, "model": "paired_contrast_ESM2_minus_ESM1", "n_variants": int(direct.loc[direct["protein"].eq(protein), "n_variants"].iloc[0]), "n_positions": int(bootstrap.loc[bootstrap["protein"].eq(protein), "sampled_positions"].iloc[0].__len__()), "prevalence": np.nan, "ap": values.mean(), "ap_ci_low": values.quantile(0.025), "ap_ci_high": values.quantile(0.975), "invalid_fraction": float(values.isna().mean()), "unstable": float(values.isna().mean()) > 0.10, "seed": seed, "n_boot": draws, "ci_interpretation": "paired measured-position bootstrap contrast", "full_draw_output": str(output_path), "full_draw_sha256": digest, "full_draw_rows": len(bootstrap)})
        rows[-1]["invalid_fraction"] = invalid_fraction
        rows[-1]["unstable"] = invalid_fraction > 0.10
    return pd.DataFrame(rows)


def _observed_permutation(eligible: pd.DataFrame, seed: int, permutations: int, output_path: Path) -> pd.DataFrame:
    """Run the observed-cohort position-block null and return compact summaries."""
    chunks = []
    for protein, group in eligible.groupby("protein", sort=True):
        null = position_block_permutation(group, seed, permutations, score_columns=MODEL_COLUMNS)
        null.insert(0, "protein", protein)
        chunks.append(null)
    draws = pd.concat(chunks, ignore_index=True)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    draws.to_csv(output_path, index=False)
    digest = _sha256(output_path)
    summaries = []
    for protein, group in draws.groupby("protein", sort=True):
        source = eligible.loc[eligible["protein"].eq(protein)]
        observed_y = source["primary_endpoint"].eq("function-impaired").astype(int).to_numpy()
        for model in MODEL_COLUMNS:
            null = group.loc[group["model"].eq(model), "ap"].dropna()
            observed_ap = evaluate_binary(observed_y, source[model].to_numpy())["ap"]
            summaries.append({"protein": protein, "model": model, "observed_ap": observed_ap, "null_mean_ap": null.mean(), "null_q025_ap": null.quantile(0.025), "null_median_ap": null.quantile(0.5), "null_q975_ap": null.quantile(0.975), "null_prevalence": group.loc[group["model"].eq(model), "prevalence"].iloc[0], "control_exceedance_fraction": float((null >= observed_ap).mean()), "seed": seed, "n_permutations": permutations, "label_block_rule": group["label_block_rule"].iloc[0], "control_interpretation": "observed-cohort null control; not a significance claim", "full_draw_output": str(output_path), "full_draw_sha256": digest, "full_draw_rows": len(draws)})
    return pd.DataFrame(summaries)


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as source:
        for chunk in iter(lambda: source.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def validate_frozen_inputs(cohort_path: Path, atlas_path: Path, provenance_path: Path) -> dict[str, str]:
    """Return input hashes and fail closed when a frozen run changes inputs."""
    input_hashes = {"cohort": _sha256(cohort_path), "atlas": _sha256(atlas_path)}
    if provenance_path.exists():
        prior = json.loads(provenance_path.read_text(encoding="utf-8"))
        if prior.get("input_sha256") and prior["input_sha256"] != input_hashes:
            raise ValueError("frozen Task 4 input SHA256 mismatch; refusing to evaluate changed inputs")
    return input_hashes


def _result_rows(eligible: pd.DataFrame, direct: pd.DataFrame, oof: list[pd.DataFrame], lopo: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for row in direct.loc[~direct["protein"].str.startswith("macro-average")].itertuples(index=False):
        values = row._asdict()
        values.update({"setting": "direct_zero_shot", "n_positions": int(eligible.loc[eligible["protein"].eq(row.protein), "position"].nunique())})
        rows.append(values)
    for prediction in oof:
        metrics = evaluate_binary(prediction["y_true"].to_numpy(), prediction["prediction"].to_numpy())
        for protein in prediction["protein"].unique():
            group = prediction.loc[prediction["protein"].eq(protein)]
            rows.append({"setting": "within_protein_nested_grouped_cv", "model": prediction["feature_set"].iloc[0], "protein": protein, **metrics, "n_positions": int(group["position"].nunique()), "evaluation": prediction.attrs.get("split", "")})
    for (protein, feature), group in lopo.groupby(["held_out_protein", "feature_set"]):
        rows.append({"setting": "leave_one_protein_out", "model": feature, "protein": protein, **evaluate_binary(group["y_true"].to_numpy(), group["prediction"].to_numpy()), "n_positions": int(group["position"].nunique()), "evaluation": lopo.attrs.get("split", "")})
    return pd.DataFrame(rows)


def _forest(direct: pd.DataFrame, bootstrap: pd.DataFrame, output: Path) -> None:
    rows = direct.loc[direct["protein"].isin(PASS_PROTEINS)].copy()
    def sampled_count(value: object) -> int:
        parsed = ast.literal_eval(value) if isinstance(value, str) else value
        return len(parsed)
    rows["n_positions"] = rows["protein"].map(
        bootstrap.groupby("protein")["sampled_positions"].first().map(sampled_count)
    )
    # Use the original eligible-row count from the direct metric table; this
    # is already the frozen common cohort denominator.
    contrasts = bootstrap.groupby("protein")["paired_ap_contrast_esm2_minus_esm1"].agg(["mean", lambda values: values.quantile(0.025), lambda values: values.quantile(0.975)]).rename(columns={"<lambda_0>": "low", "<lambda_1>": "high"})
    rows["ci_low"] = [bootstrap.loc[(bootstrap["protein"].eq(row.protein)) & (bootstrap["model"].eq(row.model)), "ap"].quantile(0.025) for row in rows.itertuples()]
    rows["ci_high"] = [bootstrap.loc[(bootstrap["protein"].eq(row.protein)) & (bootstrap["model"].eq(row.model)), "ap"].quantile(0.975) for row in rows.itertuples()]
    rows = rows.sort_values(["protein", "model"]).reset_index(drop=True)
    fig, ax = plt.subplots(figsize=(11, 5.5))
    y = np.arange(len(rows))
    for index, row in rows.iterrows():
        color = "#2c7fb8" if row["model"] == MODEL_COLUMNS[0] else "#d95f02"
        ax.errorbar(row["ap"], index, xerr=[[row["ap"] - row["ci_low"]], [row["ci_high"] - row["ap"]]], fmt="o", color=color, capsize=4, label=row["model"] if index < 2 else None)
        ax.axvline(row["prevalence"], color="#666666", linestyle="--", linewidth=0.7)
        ax.text(1.01, index, f"n={int(row['n_variants'])}, positions={int(row['n_positions'])}", va="center", transform=ax.get_yaxis_transform(), fontsize=8)
    labels = [f"{row.protein} — {row.model.replace('_impairment', '').replace('_single', '')}" for row in rows.itertuples()]
    ax.set_yticks(y, labels)
    ax.set_xlabel("Average precision (point; 95% position-bootstrap CI)")
    ax.set_title("Figure 4. Frozen zero-shot impairment benchmark\nDashed lines are per-protein prevalence baselines; CI covers measured positions only")
    ax.text(0.02, 0.98, "Blue: ESM1v | Orange: ESM2", transform=ax.transAxes, ha="left", va="top", fontsize=8)
    contrast_text = "  ".join(f"{protein}: paired ESM2−ESM1 AP contrast {values['mean']:.3f} [{values['low']:.3f}, {values['high']:.3f}]" for protein, values in contrasts.iterrows())
    fig.text(0.5, 0.01, contrast_text, ha="center", fontsize=8)
    fig.subplots_adjust(left=0.25, right=0.74, bottom=0.20, top=0.84)
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=180)
    plt.close(fig)


def _negative_control(seed: int = 2026, n: int = 10000) -> dict[str, object]:
    rng = np.random.default_rng(seed)
    y = rng.binomial(1, 0.2, n)
    random_score = rng.normal(size=n)
    result = evaluate_binary(y, random_score)
    tolerance = max(0.02, 5 / np.sqrt(n))
    result.update({"control": "large synthetic permutation/no-signal", "seed": seed, "n": n, "tolerance": tolerance, "near_prevalence": abs(result["ap"] - result["prevalence"]) <= tolerance})
    return result


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--cohort", type=Path, default=Path("D:/protein_assay_work/processed/paired_phenotype_cohort.csv"))
    parser.add_argument("--atlas", type=Path, default=external_artifact_path("assay_atlas.csv"))
    parser.add_argument("--metadata-dir", type=Path, default=Path("results/metadata"))
    parser.add_argument("--figure-dir", type=Path, default=Path("results/figures"))
    parser.add_argument("--seed", type=int, default=2026)
    parser.add_argument("--bootstrap-draws", type=int, default=2000)
    parser.add_argument("--permutation-draws", type=int, default=2000)
    args = parser.parse_args()
    args.metadata_dir.mkdir(parents=True, exist_ok=True)
    prior_provenance_path = args.metadata_dir / "task4_provenance.json"
    input_hashes = validate_frozen_inputs(args.cohort, args.atlas, prior_provenance_path)
    eligible = _read_primary(args.cohort, args.atlas)
    direct = evaluate_primary(eligible, seed=args.seed, bootstrap_draws=args.bootstrap_draws)
    bootstrap = _direct_bootstrap(eligible, args.seed, args.bootstrap_draws)
    within = [nested_grouped_cv(eligible.loc[eligible["protein"].eq(protein)].reset_index(drop=True), feature, outer_splits=5) for protein in PASS_PROTEINS for feature in ("B0", "B1", "B2")]
    within_frame = pd.concat(within, ignore_index=True)
    lopo_parts = [leave_one_protein_out(eligible, feature) for feature in ("B0", "B1", "B2")]
    lopo = pd.concat(lopo_parts, ignore_index=True)
    lopo.attrs["split"] = "leave-one-protein-out; two PASS proteins only; inner grouped-position CV within training protein"
    metrics = _result_rows(eligible, direct, within, lopo)
    metrics.to_csv(args.metadata_dir / "task4_results.csv", index=False)
    direct.to_csv(args.metadata_dir / "task4_primary_metrics.csv", index=False)
    bootstrap_summary = _bootstrap_summary(direct, bootstrap, Path("D:/protein_assay_work/processed/task4_bootstrap_draws.csv"), args.seed, args.bootstrap_draws)
    bootstrap_summary.to_csv(args.metadata_dir / "task4_bootstrap_summary.csv", index=False)
    permutation_summary = _observed_permutation(eligible, args.seed, args.permutation_draws, Path("D:/protein_assay_work/processed/task4_observed_permutation_draws.csv"))
    permutation_summary.to_csv(args.metadata_dir / "task4_observed_permutation_summary.csv", index=False)
    oof_path = external_artifact_path("task4_oof_predictions.csv")
    lopo_path = external_artifact_path("task4_lopo_predictions.csv")
    oof_path.parent.mkdir(parents=True, exist_ok=True)
    within_frame.to_csv(oof_path, index=False)
    lopo.to_csv(lopo_path, index=False)
    pd.DataFrame([_negative_control(args.seed)]).to_csv(args.metadata_dir / "task4_negative_control.csv", index=False)
    figure_path = args.figure_dir / "figure4_predictor_forest.png"
    _forest(direct, bootstrap, figure_path)
    output_paths = [args.metadata_dir / name for name in ("task4_results.csv", "task4_primary_metrics.csv", "task4_bootstrap_summary.csv", "task4_observed_permutation_summary.csv", "task4_negative_control.csv")] + [oof_path, lopo_path, figure_path]
    provenance = {"code_version": "task4-fix1", "input_sha256": input_hashes, "seed": args.seed, "bootstrap_draws": args.bootstrap_draws, "permutation_draws": args.permutation_draws, "primary_proteins": PASS_PROTEINS, "primary_endpoint": "primary_eligible=True; function-impaired versus function-preserved", "score_direction": "frozen ProteinGym impairment-oriented scores; higher means more impaired", "split_settings": {"within": "outer 5-fold grouped-position CV with 3-fold grouped-position inner selection (fallback outer 3 when required)", "leave_one_protein_out": "two PASS proteins only; inner grouped-position CV within single training protein"}, "outputs": {str(path): {"sha256": _sha256(path), "rows": int(pd.read_csv(path).shape[0]) if path.suffix == ".csv" else None} for path in output_paths}, "external_draws": {"bootstrap": {"path": "D:/protein_assay_work/processed/task4_bootstrap_draws.csv", "sha256": bootstrap_summary["full_draw_sha256"].iloc[0], "rows": len(bootstrap)}, "observed_permutation": {"path": "D:/protein_assay_work/processed/task4_observed_permutation_draws.csv", "sha256": permutation_summary["full_draw_sha256"].iloc[0], "rows": int(permutation_summary["full_draw_rows"].iloc[0])}}}
    prior_provenance_path.write_text(json.dumps(provenance, indent=2) + "\n", encoding="utf-8")
    print(metrics.to_json(orient="records", indent=2))


if __name__ == "__main__":
    main()
