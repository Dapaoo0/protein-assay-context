"""Build the paired assay atlas and audit ProteinGym zero-shot coverage."""

from __future__ import annotations

import argparse
import hashlib
import json
from datetime import UTC, datetime
from pathlib import Path
from zipfile import ZipFile

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
import yaml
from matplotlib.patches import Patch

from assay_context.artifacts import external_artifact_path
from assay_context.atlas import summarize_phenotypes
from assay_context.scores import audit_cross_assay_scores, join_scores, orient_scores

ARCHIVE_URL = "https://marks.hms.harvard.edu/proteingym/ProteinGym_v1.3/zero_shot_substitutions_scores.zip"
ARCHIVE_SIZE = 1_911_045_703
ARCHIVE_SHA256 = "3fd7cdb5e78f1d43cabfabfeb6578c252b63af23ba2ab44db0094dc3a42de36d"
MODELS = ("ESM1v_single", "ESM2_650M")
MODEL_DIRECTION = {
    # ProteinGym's zero-shot fitness scores use the native model fitness
    # orientation: higher score predicts higher fitness, so impairment is -score.
    "ESM1v_single": False,
    "ESM2_650M": False,
}
ASSAYS = {
    "CP2C9_HUMAN_Amorosi_2021_abundance.csv": "CYP2C9",
    "CP2C9_HUMAN_Amorosi_2021_activity.csv": "CYP2C9",
    "HXK4_HUMAN_Gersing_2023_abundance.csv": "GCK",
    "HXK4_HUMAN_Gersing_2022_activity.csv": "GCK",
    "PTEN_HUMAN_Matreyek_2021.csv": "PTEN",
    "PTEN_HUMAN_Mighell_2018.csv": "PTEN",
    "VKOR1_HUMAN_Chiasson_2020_abundance.csv": "VKORC1",
    "VKOR1_HUMAN_Chiasson_2020_activity.csv": "VKORC1",
}
ABUNDANCE_ASSAYS = {
    "CP2C9_HUMAN_Amorosi_2021_abundance.csv": "CYP2C9",
    "HXK4_HUMAN_Gersing_2023_abundance.csv": "GCK",
    "PTEN_HUMAN_Matreyek_2021.csv": "PTEN",
    "VKOR1_HUMAN_Chiasson_2020_abundance.csv": "VKORC1",
}
FUNCTION_ASSAYS = {
    "CP2C9_HUMAN_Amorosi_2021_activity.csv": "CYP2C9",
    "HXK4_HUMAN_Gersing_2022_activity.csv": "GCK",
    "PTEN_HUMAN_Mighell_2018.csv": "PTEN",
    "VKOR1_HUMAN_Chiasson_2020_activity.csv": "VKORC1",
}
MODEL_SOURCE_RULE = (
    "For each protein, the abundance-assay ProteinGym file is the canonical production "
    "source; the paired function-assay file is audit-only and a Task 5 sensitivity source."
)


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _write_if_changed(path: Path, text: str) -> None:
    """Avoid mutating frozen metadata on deterministic reruns."""
    encoded = text.encode("utf-8")
    if path.exists() and path.read_bytes() == encoded:
        return
    path.write_bytes(encoded)


def frozen_download_timestamp(prior: dict[str, object], archive_sha: str, now: str) -> str:
    """Keep a prior download event stable when regenerating identical inputs."""
    prior_archive = prior.get("score_archive", {})
    if isinstance(prior_archive, dict) and prior_archive.get("sha256") == archive_sha:
        timestamp = prior_archive.get("downloaded_at")
        if isinstance(timestamp, str) and timestamp:
            return timestamp
    return now


def ensure_extracted_scores(score_dir: Path, archive: Path, expected: dict[str, str]) -> dict[str, str]:
    """Validate frozen extracted hashes, or extract missing files from the verified archive."""
    score_dir.mkdir(parents=True, exist_ok=True)
    with ZipFile(archive) as zipped:
        for filename in ASSAYS:
            path = score_dir / filename
            expected_hash = expected.get(filename)
            if not path.exists():
                member = f"DMS_ProteinGym_substitutions/{filename}"
                with zipped.open(member) as source, path.open("wb") as target:
                    target.write(source.read())
            actual_hash = _sha256(path)
            if expected_hash is not None and actual_hash != expected_hash:
                raise ValueError(f"Frozen extracted score SHA256 mismatch for {filename}: {actual_hash}")
            expected[filename] = actual_hash
    return expected


def read_score_tables(
    score_dir: Path,
    metadata: dict[str, dict[str, str]],
    assays: dict[str, str] = ASSAYS,
) -> tuple[pd.DataFrame, list[dict[str, object]]]:
    """Read only the two frozen model columns from selected extracted files."""
    frames: list[pd.DataFrame] = []
    manifest: list[dict[str, object]] = []
    for filename, protein in assays.items():
        path = score_dir / filename
        if not path.exists():
            raise FileNotFoundError(f"Missing extracted score table: {path}")
        frame = pd.read_csv(path, usecols=["mutant", *MODELS])
        frame = frame.rename(columns={"mutant": "variant"})
        parsed = frame["variant"].astype("string").str.extract(r"^(?P<wt>[A-Z])(?P<position>[0-9]+)(?P<mutant>[A-Z])$")
        frame["protein"] = protein
        frame["protein_id"] = metadata[protein]["protein_id"]
        frame["sequence_hash"] = metadata[protein]["sequence_hash"]
        frame["wt"] = parsed["wt"]
        frame["position"] = pd.to_numeric(parsed["position"], errors="raise")
        frame["mutant"] = parsed["mutant"]
        frame["assay_id"] = filename.removesuffix(".csv")
        frames.append(frame.melt(id_vars=["protein", "protein_id", "sequence_hash", "variant", "wt", "position", "mutant", "assay_id"], var_name="model", value_name="score"))
        manifest.append(
            {
                "filename": filename,
                "protein": protein,
                "role": "abundance" if filename in ABUNDANCE_ASSAYS else "function",
                "rows": len(frame),
                "sha256": _sha256(path),
            }
        )
    return pd.concat(frames, ignore_index=True), manifest


def coverage_table(pairs: pd.DataFrame, joined: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for protein, group in pairs.groupby("protein", sort=True):
        covered = joined.loc[group.index]
        primary = group.get("primary_eligible", pd.Series(False, index=group.index)).astype(bool)
        for cohort_name, mask in (("all_paired", pd.Series(True, index=group.index)), ("primary_eligible", primary)):
            total = int(mask.sum())
            for model in MODELS:
                matched = int(covered.loc[mask, model].notna().sum())
                rows.append(
                    {
                        "protein": protein,
                        "cohort": cohort_name,
                        "model": model,
                        "matched": matched,
                        "total": total,
                        "coverage_fraction": matched / total if total else 0.0,
                        "gate_b": "NOT_APPLICABLE" if total == 0 else ("PASS" if matched / total >= 0.90 else "REDUCE_SCOPE"),
                    }
                )
    return pd.DataFrame(rows)


def cohort_flow_table(pairs: pd.DataFrame, joined: pd.DataFrame, attrition: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for protein, group in pairs.groupby("protein", sort=True):
        source = attrition.loc[attrition["protein"].eq(protein)].iloc[0]
        primary = group["primary_eligible"].astype(bool)
        endpoint_labeled = group.get(
            "endpoint_labeled",
            group.get("primary_endpoint", pd.Series(index=group.index, dtype="string"))
            .isin(("function-impaired", "function-preserved")),
        ).astype(bool)
        endpoint_count = int(endpoint_labeled.sum())
        primary_count = int(primary.sum())
        status = "primary" if primary_count else "descriptive-only (no primary benchmark rows)"
        stages = {
            "abundance assay input": int(source["abundance_input"]),
            "function assay input": int(source["function_input"]),
            "sequence-compatible missense": int(source["shared_proteingym_pairs"]),
            "QC-mapped paired rows": int(source["author_annotated_pairs"]),
            "endpoint labeled": endpoint_count,
            "primary benchmark eligible": primary_count,
            "ESM1v covered (primary)": int(joined.loc[primary.index[primary], "ESM1v_single"].notna().sum()),
            "both models covered (primary)": int(joined.loc[primary.index[primary], list(MODELS)].notna().all(axis=1).sum()),
        }
        rows.extend(
            {"protein": protein, "stage": stage, "count": count, "funnel": status}
            for stage, count in stages.items()
        )
    return pd.DataFrame(rows)


def assay_context_table(inventory: pd.DataFrame, assay_map: dict[str, tuple[str, str]]) -> pd.DataFrame:
    """Return compact assay IDs/readouts for the Figure 2 caption table."""
    by_id = inventory.set_index("assay_id")
    rows = []
    for protein, (abundance_id, function_id) in assay_map.items():
        rows.append(
            {
                "protein": protein,
                "abundance_assay_id": abundance_id,
                "abundance_readout": by_id.loc[abundance_id, "readout"],
                "function_assay_id": function_id,
                "function_readout": by_id.loc[function_id, "readout"],
                "numeric_threshold_note": (
                    "abundance 0.60; function loss 0.66; hyperactivity 1.18"
                    if protein == "GCK"
                    else "categorical author classes; no single numeric threshold"
                ),
            }
        )
    return pd.DataFrame(rows)


def _write_figures(
    pairs: pd.DataFrame,
    joined: pd.DataFrame,
    fractions: pd.DataFrame,
    figure_dir: Path,
    flow: pd.DataFrame,
    context: pd.DataFrame,
) -> None:
    figure_dir.mkdir(parents=True, exist_ok=True)
    stages = flow["stage"].drop_duplicates().tolist()
    stage_labels = {
        "abundance assay input": "abundance\ninput",
        "function assay input": "function\ninput",
        "sequence-compatible missense": "sequence-\ncompatible",
        "QC-mapped paired rows": "QC-mapped\npaired",
        "endpoint labeled": "endpoint\nlabeled",
        "primary benchmark eligible": "primary\nbenchmark",
        "ESM1v covered (primary)": "ESM1v\ncovered",
        "both models covered (primary)": "both models\ncovered",
    }
    proteins = sorted(flow["protein"].unique())
    fig, axes = plt.subplots(2, 2, figsize=(11, 8), squeeze=False)
    for ax, protein in zip(axes.flat, proteins, strict=True):
        data = flow.loc[flow["protein"].eq(protein)].set_index("stage").loc[stages]
        bars = ax.bar(range(len(stages)), data["count"], color=["#355c7d", "#355c7d", "#6c8ebf", "#99c1de", "#b3cde3", "#a5d6c8", "#9dc3c8", "#c9e4de"])
        ax.bar_label(bars, fmt="%d", padding=2, fontsize=7)
        ax.set_xticks(range(len(stages)), [stage_labels[stage] for stage in stages], fontsize=7)
        ax.set_ylabel("Rows")
        funnel = flow.loc[flow["protein"].eq(protein), "funnel"].iloc[0]
        ax.set_title(f"{protein} ({funnel})", fontsize=10)
    fig.suptitle("Figure 1. Per-protein cohort flow; assay inputs shown separately")
    fig.tight_layout()
    fig.savefig(figure_dir / "figure1_cohort_flow.png", dpi=180)
    plt.close(fig)

    proteins = sorted(pairs["protein"].unique())
    fig, axes = plt.subplots(2, 2, figsize=(11, 9), squeeze=False)
    phenotype_colors = {
        "preserved / impaired": "#d95f02",
        "preserved / preserved": "#1b9e77",
        "low / impaired": "#7570b3",
        "low / preserved": "#e7298a",
        "uncertain / other": "#bdbdbd",
    }
    for ax, protein in zip(axes.flat, proteins, strict=True):
        group = joined.loc[joined["protein"].eq(protein)]
        phenotype = group["abundance_label"].astype("string") + " / " + group["function_label"].astype("string")
        colors = phenotype.map(phenotype_colors).fillna(phenotype_colors["uncertain / other"])
        ax.scatter(group["DMS_score_abundance"], group["DMS_score_function"], c=colors, s=5, alpha=0.35, linewidths=0)
        ax.set_title(f"{protein} (n={len(group):,})", fontsize=10)
        ax.set_xlabel("Abundance DMS score")
        ax.set_ylabel("Function DMS score")
        if protein == "GCK":
            ax.axvline(0.6, color="#555555", linestyle="--", linewidth=0.8)
            ax.axhline(0.66, color="#555555", linestyle="--", linewidth=0.8)
            ax.axhline(1.18, color="#555555", linestyle=":", linewidth=0.8)
            ax.text(0.02, 0.04, "Numeric lines: abundance 0.60; function 0.66 / 1.18", transform=ax.transAxes, fontsize=6)
        else:
            ax.text(0.02, 0.04, "Categorical author classes; no numeric line", transform=ax.transAxes, fontsize=6)
    fig.suptitle("Figure 2. Abundance versus function (assay-specific readouts)")
    fig.legend(
        handles=[Patch(facecolor=color, label=label) for label, color in phenotype_colors.items()],
        loc="lower center", bbox_to_anchor=(0.5, 0.045), ncol=3, fontsize=7,
    )
    fig.text(
        0.5,
        0.008,
        "Colors encode abundance/function phenotype classes; full assay IDs and readouts are in figure2_assay_context.csv.",
        ha="center",
        fontsize=7,
    )
    fig.subplots_adjust(bottom=0.17, top=0.91, hspace=0.30, wspace=0.28)
    fig.savefig(figure_dir / "figure2_abundance_function_scatter.png", dpi=180)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(11, 5))
    x = range(len(fractions))
    colors = ["#4c78a8", "#f58518", "#54a24b", "#e45756", "#b0b0b0"]
    ax.bar(x, fractions["fraction"], color=[colors[PHENOTYPES.index(value)] for value in fractions["phenotype"]],
           yerr=[fractions["fraction"] - fractions["fraction_ci_low"], fractions["fraction_ci_high"] - fractions["fraction"]], capsize=3)
    labels = [f"{row.protein}\n{row.phenotype.replace('_', ' ')}\n(n={row.denominator})" for row in fractions.itertuples()]
    ax.set_xticks(list(x), labels, rotation=65, ha="right", fontsize=7)
    ax.set_ylabel("Fraction of all paired rows")
    ax.set_ylim(0, 1)
    ax.set_title("Figure 3. Phenotype fractions with 95% position-bootstrap CIs")
    fig.tight_layout()
    fig.savefig(figure_dir / "figure3_phenotype_fractions.png", dpi=180)
    plt.close(fig)


PHENOTYPES = (
    "abundance-preserved_function-impaired",
    "abundance-preserved_function-preserved",
    "abundance-low_function-impaired",
    "abundance-low_function-preserved",
    "uncertain-or-other",
)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--cohort", type=Path, default=Path("D:/protein_assay_work/processed/paired_phenotype_cohort.csv"))
    parser.add_argument("--score-dir", type=Path, default=Path("D:/protein_assay_work/processed/model_scores"))
    parser.add_argument("--archive", type=Path, default=Path("D:/protein_assay_work/downloads/zero_shot_substitutions_scores.zip"))
    parser.add_argument("--metadata-dir", type=Path, default=Path("results/metadata"))
    parser.add_argument("--figure-dir", type=Path, default=Path("results/figures"))
    parser.add_argument("--config", type=Path, default=Path("config/analysis_freeze.yaml"))
    args = parser.parse_args()
    if args.archive.stat().st_size != ARCHIVE_SIZE:
        raise ValueError(f"Score archive has unexpected size: {args.archive.stat().st_size}")
    archive_sha = _sha256(args.archive)
    if archive_sha != ARCHIVE_SHA256:
        raise ValueError(f"Score archive SHA256 mismatch: {archive_sha}")
    provenance_path = args.metadata_dir / "task2_provenance.json"
    prior_provenance = json.loads(provenance_path.read_text(encoding="utf-8"))
    frozen_files = {
        row["filename"]: row["sha256"]
        for row in prior_provenance.get("score_archive", {}).get("extracted_target_files", [])
        if "filename" in row and "sha256" in row
    }
    ensure_extracted_scores(args.score_dir, args.archive, frozen_files)
    inventory = pd.read_csv("results/metadata/assay_inventory.csv")
    metadata_by_protein = {
        protein_id: {"protein_id": protein_id, "sequence_hash": str(group.iloc[0]["sequence_hash"])}
        for protein_id, group in inventory.groupby("protein_id", sort=False)
    }
    protein_metadata = {
        "CYP2C9": metadata_by_protein["CP2C9_HUMAN"],
        "GCK": metadata_by_protein["HXK4_HUMAN"],
        "PTEN": metadata_by_protein["PTEN_HUMAN"],
        "VKORC1": metadata_by_protein["VKOR1_HUMAN"],
    }
    pairs = pd.read_csv(args.cohort)
    pairs["protein_id"] = pairs["protein"].map({protein: value["protein_id"] for protein, value in protein_metadata.items()})
    expected_hashes = pairs["protein"].map({protein: value["sequence_hash"] for protein, value in protein_metadata.items()})
    if pairs["protein_id"].isna().any() or ~pairs["sequence_hash"].eq(expected_hashes).all():
        raise ValueError("Cohort protein_id/sequence_hash does not match assay inventory")
    assay_map = {
        "CYP2C9": ("CP2C9_HUMAN_Amorosi_2021_abundance", "CP2C9_HUMAN_Amorosi_2021_activity"),
        "PTEN": ("PTEN_HUMAN_Matreyek_2021", "PTEN_HUMAN_Mighell_2018"),
        "GCK": ("HXK4_HUMAN_Gersing_2023_abundance", "HXK4_HUMAN_Gersing_2022_activity"),
        "VKORC1": ("VKOR1_HUMAN_Chiasson_2020_abundance", "VKOR1_HUMAN_Chiasson_2020_activity"),
    }
    pairs["abundance_assay_id"] = pairs["protein"].map({key: value[0] for key, value in assay_map.items()})
    pairs["function_assay_id"] = pairs["protein"].map({key: value[1] for key, value in assay_map.items()})
    abundance_scores, abundance_manifest = read_score_tables(
        args.score_dir, protein_metadata, ABUNDANCE_ASSAYS
    )
    function_scores, function_manifest = read_score_tables(
        args.score_dir, protein_metadata, FUNCTION_ASSAYS
    )
    # Production scores come from one uniformly selected abundance file per
    # protein. Function scores are never coalesced into this production table.
    joined = join_scores(pairs, abundance_scores)
    joined = orient_scores(joined, MODEL_DIRECTION)
    function_joined = orient_scores(join_scores(pairs, function_scores), MODEL_DIRECTION)
    coverage = coverage_table(pairs, joined)
    primary_coverage = coverage.loc[
        coverage["cohort"].eq("primary_eligible") & coverage["total"].gt(0)
    ]
    if not primary_coverage["gate_b"].eq("PASS").all():
        primary_contrast = "REDUCE_SCOPE: one or more Gate A proteins lacks >=90% coverage"
    else:
        primary_contrast = "ESM1v_single_vs_ESM2_650M"
    fractions = summarize_phenotypes(pairs)
    attrition = pd.read_csv(args.metadata_dir / "task2_attrition.csv")
    flow = cohort_flow_table(pairs, joined, attrition)
    context = assay_context_table(inventory, assay_map)
    args.metadata_dir.mkdir(parents=True, exist_ok=True)
    coverage.to_csv(args.metadata_dir / "score_coverage.csv", index=False)
    function_coverage = coverage_table(pairs, function_joined)
    function_coverage.insert(1, "source", "function_sensitivity")
    function_coverage.to_csv(args.metadata_dir / "function_score_coverage.csv", index=False)
    detail_path = Path("D:/protein_assay_work/processed/score_consistency_detail.csv")
    detail_path.parent.mkdir(parents=True, exist_ok=True)
    summary_rows: list[dict[str, object]] = []
    detail_frames: list[pd.DataFrame] = []
    for protein in sorted(pairs["protein"].unique()):
        for model in MODELS:
            summary, model_detail = audit_cross_assay_scores(
                abundance_scores.loc[abundance_scores["protein"].eq(protein) & abundance_scores["model"].eq(model)],
                function_scores.loc[function_scores["protein"].eq(protein) & function_scores["model"].eq(model)],
            )
            summary_rows.append({"protein": protein, "model": model, **summary})
            detail_frames.append(model_detail)
    detail_columns = [
        "protein_id", "sequence_hash", "wt", "position", "mutant", "protein", "variant", "model",
        "assay_id_1", "value_1", "assay_id_2", "value_2", "absolute_difference",
    ]
    detail = pd.concat(detail_frames, ignore_index=True) if detail_frames else pd.DataFrame(columns=detail_columns)
    detail = detail.reindex(columns=detail_columns)
    detail.to_csv(detail_path, index=False)
    pd.DataFrame(summary_rows).to_csv(args.metadata_dir / "score_consistency.csv", index=False)
    flow.to_csv(args.metadata_dir / "cohort_flow.csv", index=False)
    context.to_csv(args.metadata_dir / "figure2_assay_context.csv", index=False)
    fractions.to_csv(args.metadata_dir / "phenotype_fractions.csv", index=False)
    atlas_columns = ["protein", "protein_id", "sequence_hash", "variant", "wt", "position", "mutant", "abundance_assay_id", "function_assay_id", "abundance_label", "function_label", "primary_endpoint", "endpoint_labeled", *MODELS, *[f"{m}_impairment" for m in MODELS]]
    atlas_path = external_artifact_path("assay_atlas.csv")
    atlas_path.parent.mkdir(parents=True, exist_ok=True)
    joined.loc[:, [column for column in atlas_columns if column in joined]].to_csv(atlas_path, index=False)
    _write_figures(pairs, joined, fractions, args.figure_dir, flow, context)
    function_sensitivity_path = Path("D:/protein_assay_work/processed/function_sensitivity_scores.csv")
    function_sensitivity_path.parent.mkdir(parents=True, exist_ok=True)
    function_scores.to_csv(function_sensitivity_path, index=False)
    timestamp = frozen_download_timestamp(prior_provenance, archive_sha, datetime.now(UTC).isoformat())
    archive_manifest = {
        "url": ARCHIVE_URL,
        "path": str(args.archive),
        "content_length": ARCHIVE_SIZE,
        "size_bytes": args.archive.stat().st_size,
        "sha256": archive_sha,
        "downloaded_at": timestamp,
        "extracted_target_files": abundance_manifest + function_manifest,
        "cross_assay_discrepancy_detail": {
            "path": str(detail_path),
            "rows": len(detail),
            "sha256": _sha256(detail_path),
            "rule": "exact numeric comparison; non-identical values are upstream cross-assay numerical discrepancies",
        },
        "function_sensitivity_source": {
            "path": str(function_sensitivity_path),
            "rows": len(function_scores),
            "sha256": _sha256(function_sensitivity_path),
        },
    }
    freeze = yaml.safe_load(args.config.read_text(encoding="utf-8"))
    freeze["score_archive"] = archive_manifest
    freeze["frozen_models"] = {
        model: {"column": model, "higher_is_impaired": MODEL_DIRECTION[model], "source": "ProteinGym v1.3 zero-shot score column"}
        for model in MODELS
    }
    freeze["model_evaluation"] = {
        "gate_b": "PASS" if primary_coverage["gate_b"].eq("PASS").all() else "REDUCE_SCOPE",
        "primary_contrast": primary_contrast,
        "coverage_table": "results/metadata/score_coverage.csv",
        "direction_rule": "ProteinGym zero-shot score columns are native fitness orientation; impairment score is negative raw score",
        "production_score_source": MODEL_SOURCE_RULE,
        "cross_assay_audit_rule": "exact numeric comparison; non-identical values are reported as upstream cross-assay numerical discrepancies and are not coalesced",
        "function_sensitivity_source": "D:\\protein_assay_work\\processed\\function_sensitivity_scores.csv",
    }
    _write_if_changed(args.config, yaml.safe_dump(freeze, sort_keys=False, allow_unicode=False))
    provenance_path = args.metadata_dir / "task2_provenance.json"
    provenance = prior_provenance
    provenance["score_archive"] = archive_manifest
    provenance["frozen_models"] = freeze["frozen_models"]
    provenance["score_coverage"] = coverage.to_dict(orient="records")
    provenance["primary_model_contrast"] = primary_contrast
    _write_if_changed(provenance_path, json.dumps(provenance, indent=2) + "\n")
    print(json.dumps({"rows": len(joined), "coverage": coverage.to_dict(orient="records"), "primary_contrast": primary_contrast}, indent=2))


if __name__ == "__main__":
    main()
