"""Create the paired-assay inventory and mapping feasibility audit."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import pandas as pd
import yaml

from assay_context.inventory import audit_pair_overlap, build_inventory, validate_pair_config


def file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--reference", type=Path, required=True)
    parser.add_argument("--archive", type=Path, required=True)
    parser.add_argument("--assay-dir", type=Path, required=True)
    parser.add_argument("--config", type=Path, default=Path("config/candidates.yaml"))
    parser.add_argument("--output-dir", type=Path, default=Path("results/metadata"))
    return parser.parse_args()


def _mapping_gate(audit: dict[str, int | bool]) -> tuple[str, str]:
    """Return a machine-readable mapping status and reason."""
    checks = (
        (audit["shared_missense"] >= 300, "shared canonical missense variants < 300"),
        (audit["shared_positions"] >= 50, "shared canonical positions < 50"),
        (audit["bad_wild_type"] == 0, "wild-type/reference mismatches detected"),
    )
    failures = [reason for passed, reason in checks if not passed]
    if failures:
        return "REDUCE_SCOPE", "; ".join(failures)
    return "PASS", "Exact canonical missense overlap meets mapping thresholds"


def _require_assay_file(path: Path, pair_name: str, role: str) -> None:
    if not path.is_file():
        raise ValueError(f"Pair {pair_name} {role} assay file not found: {path}")


def overlap_output_row(
    pair: dict[str, str],
    audit: dict[str, int | bool],
    *,
    mapping_gate: str,
    mapping_gate_reason: str,
    abundance_file_sha256: str,
    function_file_sha256: str,
) -> dict[str, object]:
    """Build one overlap record using the configured cohort-facing pair name."""
    return {
        "protein": pair["pair_name"],
        "abundance_assay": pair["abundance_assay"],
        "function_assay": pair["function_assay"],
        **audit,
        "mapping_gate": mapping_gate,
        "mapping_gate_reason": mapping_gate_reason,
        "abundance_file_sha256": abundance_file_sha256,
        "function_file_sha256": function_file_sha256,
    }


def unique_reasons(reasons: list[str]) -> list[str]:
    """Deduplicate repeated gate reasons without changing their order."""
    return list(dict.fromkeys(reasons))


def main() -> None:
    args = parse_args()
    config = yaml.safe_load(args.config.read_text(encoding="utf-8"))
    source = config["source"]
    reference_hash = file_sha256(args.reference)
    archive_hash = file_sha256(args.archive)
    if reference_hash != source["sha256"]:
        raise ValueError(f"Reference SHA256 mismatch: {reference_hash}")
    if archive_hash != source["archive_sha256"]:
        raise ValueError(f"Archive SHA256 mismatch: {archive_hash}")

    reference = pd.read_csv(args.reference)
    inventory = build_inventory(reference, config)
    pair_rows = validate_pair_config(inventory, config.get("pairs"))
    by_id = inventory.set_index("assay_id", verify_integrity=True)
    overlaps: list[dict[str, object]] = []
    for pair in pair_rows:
        pair_name = pair["pair_name"]
        abundance_id = pair["abundance_assay"]
        function_id = pair["function_assay"]
        abundance_path = args.assay_dir / by_id.at[abundance_id, "reference_filename"]
        function_path = args.assay_dir / by_id.at[function_id, "reference_filename"]
        _require_assay_file(abundance_path, pair_name, "abundance")
        _require_assay_file(function_path, pair_name, "function")
        abundance = pd.read_csv(abundance_path)
        function = pd.read_csv(function_path)
        audit = audit_pair_overlap(
            abundance,
            function,
            by_id.at[abundance_id, "sequence"],
            by_id.at[function_id, "sequence"],
        )
        mapping_gate, mapping_gate_reason = _mapping_gate(audit)
        overlaps.append(
            overlap_output_row(
                pair,
                audit,
                mapping_gate=mapping_gate,
                mapping_gate_reason=mapping_gate_reason,
                abundance_file_sha256=file_sha256(abundance_path),
                function_file_sha256=file_sha256(function_path),
            )
        )

    args.output_dir.mkdir(parents=True, exist_ok=True)
    inventory.to_csv(args.output_dir / "assay_inventory.csv", index=False)
    overlap_frame = pd.DataFrame(overlaps)
    overlap_frame.to_csv(args.output_dir / "pair_overlap.csv", index=False)
    mapping_status = (
        "PASS" if overlap_frame["mapping_gate"].eq("PASS").all() else "REDUCE_SCOPE"
    )
    mapping_reasons = unique_reasons(overlap_frame["mapping_gate_reason"].tolist())
    class_balance_reason = (
        "REDUCE_SCOPE until Task 2 freezes author-calibrated phenotype thresholds, "
        "uncertainty handling, and the minimum class/position counts"
    )
    model_score_reason = (
        "REDUCE_SCOPE until Task 3 audits pretrained-score coverage and freezes the score cohort"
    )
    overall_reasons = list(mapping_reasons) + [class_balance_reason, model_score_reason]
    summary = {
        "reference_sha256": reference_hash,
        "archive_sha256": archive_hash,
        "proteins_audited": len(overlaps),
        "mapping_gate_passes": int(overlap_frame["mapping_gate"].eq("PASS").sum()),
        "mapping_gate": mapping_status,
        "mapping_gate_reasons": mapping_reasons,
        "class_balance_gate": "REDUCE_SCOPE",
        "class_balance_reasons": [class_balance_reason],
        "model_score_gate": "REDUCE_SCOPE",
        "model_score_reasons": [model_score_reason],
        "overall_verdict": "REDUCE_SCOPE",
        "overall_reasons": overall_reasons,
        "gates": {
            "mapping": {"status": mapping_status, "reasons": mapping_reasons},
            "class_balance": {"status": "REDUCE_SCOPE", "reasons": [class_balance_reason]},
            "model_score": {"status": "REDUCE_SCOPE", "reasons": [model_score_reason]},
        },
    }
    (args.output_dir / "feasibility_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8"
    )


if __name__ == "__main__":
    main()
