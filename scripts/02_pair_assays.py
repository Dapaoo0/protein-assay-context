"""Freeze validated author-calibrated phenotype pairs for downstream analysis."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
from typing import Any

import pandas as pd
import yaml

from assay_context.phenotypes import (
    adapt_cyp2c9,
    adapt_gck,
    adapt_pten_abundance,
    adapt_pten_function,
    adapt_vkor,
    gate_a_status,
    primary_binary_eligibility,
)
from assay_context.variants import pair_assays, parse_variant

SEED = 2026
AUTHOR_SHA256 = {
    "CYP2C9_activity_abundance_scores.csv": "283883889947e47c47b08f2f096feaf7ad5dd49d0f79432b8675c10150b61124",
    "PTEN_Composite_abundance_data.tsv": "e8172400de12c53043556f3d0b07e44e5a23c21582b14f3bc8f6d3dad9c0fe6c",
    "PTEN_Mighell_phosphatase.csv": "6d0bc5e0f6a50f59ca272a19b5756107a089473d2d8e54a706e204da27601ecd",
    "GCK_abundance.csv": "75fdeece9cc1c1a35a417956b1e9d5c62c39fe0a4a7c91a364bc43111c837dc7",
    "GCK_activity.csv": "f7efb1d51dd380fcd275c56da70592634b945579bbdf3acc065ffe226c19e7cd",
    "VKOR_combined.csv": "ea4bf42a40f1a7965667eca37318b2f17fc34ae8174f656939d5626e56523a8b",
}
AUTHOR_REVISIONS = {
    "CYP2C9_activity_abundance_scores.csv": "dunhamlab/CYP2C9@e727ffe1e4746311d44c09b4f6ffa35ca3acf665",
    "PTEN_Composite_abundance_data.tsv": "matreyeklab/pten_composite@73443cfa612ae6e570347f96461ed71a2641608a",
    "PTEN_Mighell_phosphatase.csv": "matreyeklab/pten_composite@73443cfa612ae6e570347f96461ed71a2641608a",
    "GCK_abundance.csv": "KULL-Centre/_2024_Gersing_GCKabundance@e77eeb3a3a3ec89baf4ac0e4196a041f7d5d32ae",
    "GCK_activity.csv": "KULL-Centre/_2024_Gersing_GCKabundance@e77eeb3a3a3ec89baf4ac0e4196a041f7d5d32ae",
    "VKOR_combined.csv": "FowlerLab/VKOR@f80bb91695967393adacd7b33ded6632283d42f1",
}


def file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _canonical_assay(frame: pd.DataFrame, sequence: str, role: str) -> tuple[pd.DataFrame, int]:
    """Validate assay mutations against its sequence and keep missense rows."""
    if "mutant" not in frame:
        raise ValueError(f"{role} ProteinGym table is missing mutant column")
    rows: list[dict[str, Any]] = []
    excluded = 0
    for _, row in frame.iterrows():
        value = row["mutant"]
        try:
            wt, position, mutant = parse_variant(value, sequence)
        except (ValueError, TypeError):
            excluded += 1
            continue
        item = row.to_dict()
        item.update({"variant": str(value).strip().upper(), "wt": wt, "position": position, "mutant": mutant})
        rows.append(item)
    result = pd.DataFrame(rows)
    if result.empty:
        result = pd.DataFrame(columns=["variant", "wt", "position", "mutant"])
    if result.duplicated(["wt", "position", "mutant"]).any():
        raise ValueError(f"{role} ProteinGym table contains duplicate canonical variant keys")
    return result, excluded


def _author_variants(frame: pd.DataFrame, sequence: str, role: str) -> tuple[pd.DataFrame, int]:
    rows: list[dict[str, Any]] = []
    excluded = 0
    for _, row in frame.iterrows():
        try:
            wt, position, mutant = parse_variant(row["variant"], sequence)
        except (ValueError, TypeError):
            excluded += 1
            continue
        item = row.to_dict()
        item.update({"wt": wt, "position": position, "mutant": mutant})
        rows.append(item)
    result = pd.DataFrame(rows)
    if result.empty:
        result = pd.DataFrame(columns=["variant", "wt", "position", "mutant"])
    if result.duplicated(["wt", "position", "mutant"]).any():
        raise ValueError(f"{role} author table contains duplicate canonical variant keys")
    return result, excluded


def _read_author(author_dir: Path, protein: str) -> tuple[pd.DataFrame, pd.DataFrame | None, list[str]]:
    if protein == "CYP2C9":
        path = author_dir / "CYP2C9_activity_abundance_scores.csv"
        return pd.read_csv(path), None, [path.name]
    if protein == "PTEN":
        abundance = author_dir / "PTEN_Composite_abundance_data.tsv"
        function = author_dir / "PTEN_Mighell_phosphatase.csv"
        return pd.read_csv(abundance, sep="\t"), pd.read_csv(function), [abundance.name, function.name]
    if protein == "GCK":
        abundance = author_dir / "GCK_abundance.csv"
        function = author_dir / "GCK_activity.csv"
        return pd.read_csv(abundance), pd.read_csv(function), [abundance.name, function.name]
    if protein == "VKORC1":
        path = author_dir / "VKOR_combined.csv"
        return pd.read_csv(path), None, [path.name]
    raise ValueError(f"No author adapter configured for {protein}")


def _raw_author_stats(frame: pd.DataFrame, sequence: str) -> tuple[int, int, int]:
    """Count source rows before adapter filtering, including exclusions."""
    valid = 0
    for value in frame[_variant_column(frame)].tolist():
        try:
            parse_variant(value, sequence)
        except (ValueError, TypeError):
            continue
        valid += 1
    return len(frame), valid, len(frame) - valid


def _variant_column(frame: pd.DataFrame) -> str:
    for column in ("variant", "Variant (one letter)"):
        if column in frame:
            return column
    raise ValueError("Author source table is missing a variant column")


def _merge_author_metadata(abundance: pd.DataFrame, function: pd.DataFrame) -> pd.DataFrame:
    abundance = abundance.rename(columns={"uncertainty_kind": "abundance_uncertainty_kind", "source_adapter": "abundance_source_adapter"})
    function = function.rename(columns={"uncertainty_kind": "function_uncertainty_kind", "source_adapter": "function_source_adapter"})
    return abundance.merge(function, on="variant", how="outer", validate="one_to_one")


def _build_audit_sample(
    protein: str,
    raw_frames: list[tuple[str, pd.DataFrame]],
    author_adapted: pd.DataFrame,
    cohort: pd.DataFrame,
    sequence: str,
) -> pd.DataFrame:
    """Compare a fixed sample to source fields, retaining excluded rows."""
    records: list[dict[str, Any]] = []
    adapted = author_adapted.set_index("variant", verify_integrity=True)
    cohort_by_variant = {
        f"{row.wt}{row.position}{row.mutant}": row
        for row in cohort.itertuples()
    }
    for source_name, frame in raw_frames:
        variant_column = _variant_column(frame)
        for _, row in frame.iterrows():
            variant = str(row[variant_column]).strip().upper()
            try:
                parse_variant(variant, sequence)
                source_status = "mapped_missense" if variant in adapted.index else "excluded_by_adapter"
            except (ValueError, TypeError):
                source_status = "excluded_invalid_or_non_missense"
            expected = adapted.loc[variant] if variant in adapted.index else None
            cohort_row = cohort_by_variant.get(variant)
            source_abundance = row.get("abundance_class", pd.NA)
            source_function = row.get("activity_class", row.get("Cum_score", row.get("activity_score", pd.NA)))
            records.append({
                "protein": protein,
                "source_table": source_name,
                "variant": variant,
                "source_class": row.get("class", row.get("Type", pd.NA)),
                "source_abundance_field": source_abundance,
                "source_function_field": source_function,
                "source_row_status": source_status,
                "expected_abundance_label": expected.get("abundance_label", pd.NA) if expected is not None else pd.NA,
                "expected_function_label": expected.get("function_label", pd.NA) if expected is not None else pd.NA,
                "cohort_abundance_label": getattr(cohort_row, "abundance_label", pd.NA) if cohort_row is not None else pd.NA,
                "cohort_function_label": getattr(cohort_row, "function_label", pd.NA) if cohort_row is not None else pd.NA,
                "audit_seed": SEED,
            })
    source_audit = pd.DataFrame.from_records(records)
    sample = source_audit.sample(n=min(20, len(source_audit)), random_state=SEED).copy()
    sample["labels_match"] = (
        sample["source_row_status"].eq("mapped_missense")
        & sample["cohort_abundance_label"].notna()
        & sample["cohort_function_label"].notna()
        & sample["expected_abundance_label"].astype("string").fillna("__NA__").eq(
            sample["cohort_abundance_label"].astype("string").fillna("__NA__")
        )
        & sample["expected_function_label"].astype("string").fillna("__NA__").eq(
            sample["cohort_function_label"].astype("string").fillna("__NA__")
        )
    )
    sample["audit_outcome"] = "excluded"
    sample.loc[sample["source_row_status"].eq("mapped_missense") & sample["labels_match"], "audit_outcome"] = "match"
    sample.loc[sample["source_row_status"].eq("mapped_missense") & ~sample["labels_match"], "audit_outcome"] = "mismatch_or_not_in_cohort"
    return sample


def build_cohort(config: dict[str, Any], inventory: pd.DataFrame, assay_dir: Path, author_dir: Path) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame, dict[str, Any]]:
    by_id = inventory.set_index("assay_id", verify_integrity=True)
    cohorts: list[pd.DataFrame] = []
    attrition: list[dict[str, Any]] = []
    counts: list[dict[str, Any]] = []
    samples: list[pd.DataFrame] = []
    manifest: dict[str, Any] = {"author_sources": {}, "assay_sources": {}}
    for protein, pair in config["pairs"].items():
        abundance_id, function_id = pair["abundance"], pair["function"]
        abundance_meta = by_id.loc[abundance_id]
        function_meta = by_id.loc[function_id]
        if abundance_meta["protein_id"] != function_meta["protein_id"]:
            raise ValueError(f"{protein} assay protein IDs differ")
        for role, metadata in (("abundance", abundance_meta), ("function", function_meta)):
            sequence_hash = hashlib.sha256(str(metadata["sequence"]).encode("ascii")).hexdigest()
            if sequence_hash != metadata["sequence_hash"]:
                raise ValueError(f"{protein} {role} metadata sequence hash is invalid")
        if abundance_meta["sequence"] != function_meta["sequence"] or abundance_meta["sequence_hash"] != function_meta["sequence_hash"]:
            raise ValueError(f"{protein} abundance/function metadata sequences differ")
        abundance_sequence = abundance_meta["sequence"]
        function_sequence = function_meta["sequence"]
        abundance_path = assay_dir / abundance_meta["reference_filename"]
        function_path = assay_dir / function_meta["reference_filename"]
        abundance_raw = pd.read_csv(abundance_path)
        function_raw = pd.read_csv(function_path)
        abundance_c, abundance_excluded = _canonical_assay(abundance_raw, abundance_sequence, f"{protein} abundance")
        function_c, function_excluded = _canonical_assay(function_raw, function_sequence, f"{protein} function")
        pg_pairs = pair_assays(abundance_c, function_c)
        first_author, second_author, author_names = _read_author(author_dir, protein)
        if protein == "CYP2C9":
            author_adapted = adapt_cyp2c9(first_author)
            raw_author_frames = [(author_names[0], first_author)]
        elif protein == "PTEN":
            assert second_author is not None
            author_adapted = _merge_author_metadata(adapt_pten_abundance(first_author), adapt_pten_function(second_author))
            raw_author_frames = [(author_names[0], first_author), (author_names[1], second_author)]
        elif protein == "GCK":
            assert second_author is not None
            author_adapted = _merge_author_metadata(adapt_gck(first_author, role="abundance"), adapt_gck(second_author, role="function"))
            raw_author_frames = [(author_names[0], first_author), (author_names[1], second_author)]
        else:
            author_adapted = adapt_vkor(first_author)
            raw_author_frames = [(author_names[0], first_author)]
        author_c, author_excluded = _author_variants(author_adapted, abundance_sequence, f"{protein}")
        merged = pg_pairs.merge(author_c, on=["wt", "position", "mutant"], how="left", validate="one_to_one", suffixes=("", "_author"))
        # PTEN/GCK have one metadata flag per source adapter.  Keep the
        # canonical cohort-level flag and avoid ambiguous pandas suffixes in
        # the frozen artifact.
        merged.insert(0, "protein", protein)
        merged["sequence_hash"] = abundance_meta["sequence_hash"]
        merged["primary_endpoint"] = "not_eligible"
        eligible = merged["abundance_label"].eq("preserved")
        merged.loc[eligible & merged["function_label"].eq("impaired"), "primary_endpoint"] = "function-impaired"
        merged.loc[eligible & merged["function_label"].eq("preserved"), "primary_endpoint"] = "function-preserved"
        merged["endpoint_labeled"] = merged["primary_endpoint"].isin(["function-impaired", "function-preserved"])
        pair_gate = gate_a_status(merged)
        merged["gate_a"] = pair_gate
        merged["primary_eligible"] = primary_binary_eligibility(merged, pair_gate)
        cohorts.append(merged)
        author_stats = [_raw_author_stats(frame, abundance_sequence) for _, frame in raw_author_frames]
        attrition.append({
            "protein": protein,
            "abundance_input": len(abundance_raw),
            "function_input": len(function_raw),
            "abundance_invalid_or_non_missense": abundance_excluded,
            "function_invalid_or_non_missense": function_excluded,
            "abundance_valid_missense": len(abundance_c),
            "function_valid_missense": len(function_c),
            "shared_proteingym_pairs": len(pg_pairs),
            "author_source_input": int(sum(stats[0] for stats in author_stats)),
            "author_source_valid_missense": int(sum(stats[1] for stats in author_stats)),
            "author_invalid_or_non_missense": int(sum(stats[2] for stats in author_stats)),
            "author_adapter_excluded": author_excluded,
            "author_annotated_pairs": int((merged["abundance_label"].notna() & merged["function_label"].notna()).sum()),
            "primary_eligible": int(merged["primary_eligible"].sum()),
        })
        endpoint_frame = merged.loc[merged["endpoint_labeled"]]
        counts.append({
            "protein": protein,
            "paired_missense": len(merged),
            "abundance_preserved": int(merged["abundance_label"].eq("preserved").sum()),
            "function_impaired_in_abundance_preserved": int((merged["primary_endpoint"] == "function-impaired").sum()),
            "function_preserved_in_abundance_preserved": int((merged["primary_endpoint"] == "function-preserved").sum()),
            "function_uncertain_in_abundance_preserved": int((eligible & merged["function_label"].eq("uncertain")).sum()),
            "function_other_in_abundance_preserved": int((eligible & merged["function_label"].eq("other")).sum()),
            "impaired_positions": int(endpoint_frame.loc[endpoint_frame["primary_endpoint"] == "function-impaired", "position"].nunique()),
            "preserved_positions": int(endpoint_frame.loc[endpoint_frame["primary_endpoint"] == "function-preserved", "position"].nunique()),
            "gate_a": pair_gate,
        })
        samples.append(_build_audit_sample(protein, raw_author_frames, author_adapted, merged, abundance_sequence))
        for name in (abundance_path.name, function_path.name):
            manifest["assay_sources"][name] = {"path": str(assay_dir / name), "sha256": file_sha256(assay_dir / name)}
        for name in author_names:
            path = author_dir / name
            manifest["author_sources"][name] = {
                "path": str(path),
                "sha256": file_sha256(path),
                "expected_sha256": AUTHOR_SHA256[name],
                "revision": AUTHOR_REVISIONS[name],
            }
            if manifest["author_sources"][name]["sha256"] != AUTHOR_SHA256[name]:
                raise ValueError(f"Author source SHA256 mismatch for {name}")
    return pd.concat(cohorts, ignore_index=True), pd.DataFrame(attrition), pd.DataFrame(counts), pd.concat(samples, ignore_index=True), {**manifest, "audit_seed": SEED}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", type=Path, default=Path("config/candidates.yaml"))
    parser.add_argument("--inventory", type=Path, default=Path("results/metadata/assay_inventory.csv"))
    parser.add_argument("--assay-dir", type=Path, default=Path("D:/protein_assay_work/assays"))
    parser.add_argument("--author-dir", type=Path, default=Path("D:/protein_assay_work/author_source"))
    parser.add_argument("--processed-dir", type=Path, default=Path("D:/protein_assay_work/processed"))
    parser.add_argument("--output-dir", type=Path, default=Path("results/metadata"))
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    config = yaml.safe_load(args.config.read_text(encoding="utf-8"))
    inventory = pd.read_csv(args.inventory)
    cohort, attrition, counts, audit, manifest = build_cohort(config, inventory, args.assay_dir, args.author_dir)
    args.processed_dir.mkdir(parents=True, exist_ok=True)
    cohort_path = args.processed_dir / "paired_phenotype_cohort.csv"
    cohort.to_csv(cohort_path, index=False)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    attrition.to_csv(args.output_dir / "task2_attrition.csv", index=False)
    counts.to_csv(args.output_dir / "task2_counts.csv", index=False)
    audit.to_csv(args.output_dir / "task2_audit_sample.csv", index=False)
    manifest["processed_cohort"] = {"path": str(cohort_path), "sha256": file_sha256(cohort_path), "rows": len(cohort)}
    (args.output_dir / "task2_provenance.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    freeze = {
        "freeze_date": "2026-09-05",
        "task": "Task 2 paired phenotype cohort",
        "primary_endpoint": "function-impaired versus function-preserved among abundance-preserved variants",
        "proteins": sorted(counts["protein"].tolist()),
        "assays": config["pairs"],
        "author_sources": AUTHOR_SHA256,
        "thresholds": {
            "CYP2C9": {"abundance_preserved": "abundance_class == wt-like", "abundance_low": "abundance_class in {decreased,nonsense-like}", "function_impaired": "activity_class in {decreased,nonsense-like}", "uncertainty": "possibly_*, increased, or missing; author point-estimate class"},
            "PTEN": {"abundance_preserved": "abundance_class == wt-like", "function_impaired": "High_conf == True and Cum_score <= -1.11", "function_preserved": "High_conf == True and -1.11 < Cum_score <= 0.89"},
            "GCK": {"confidence_interval": "score +/- t(0.975, df) * SE", "abundance_threshold": 0.6, "function_loss_threshold": 0.66, "function_hyperactivity_threshold": 1.18, "equality": "uncertain"},
            "VKORC1": {"abundance_preserved": "abundance_class == wt-like", "function_impaired": "activity_class == low", "uncertainty": "possible/high/missing are not preserved"},
        },
        "excluded": "non-missense and invalid reference rows are retained in attrition counts, not cohort",
        "audit_seed": SEED,
        "gate_a_minimums": {"paired_missense": 300, "total_positions": 50, "variants_per_class": 50, "positions_per_class": 10},
        "primary_benchmark_requires": "gate_a == PASS and endpoint_labeled == true",
        "model_evaluation": "deferred; no model scores included in this freeze",
    }
    args_config = args.config.parent / "analysis_freeze.yaml"
    args_config.write_text(yaml.safe_dump(freeze, sort_keys=False, allow_unicode=False), encoding="utf-8")


if __name__ == "__main__":
    main()
