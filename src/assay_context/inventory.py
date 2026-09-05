"""Build a sequence-audited inventory of preselected protein assays."""

from __future__ import annotations

import hashlib
import json
import re
from collections.abc import Mapping
from typing import Any

import pandas as pd

REFERENCE_COLUMNS = (
    "DMS_id",
    "DMS_filename",
    "DMS_number_single_mutants",
    "UniProt_ID",
    "target_seq",
    "jo",
    "taxon",
    "source_organism",
    "selection_assay",
    "selection_type",
    "raw_DMS_directionality",
    "DMS_binarization_method",
    "DMS_binarization_cutoff",
)

CANDIDATE_COLUMNS = (
    "role",
    "assay_url",
    "construct",
    "host_system",
    "tags",
    "normalization",
    "controls",
    "qc_rule",
    "threshold_source",
    "license",
)

SOURCE_COLUMNS = ("name", "url", "downloaded_at", "sha256", "revision")
MISSENSE_PATTERN = re.compile(r"^([A-Z])(\d+)([A-Z])$")
STANDARD_AMINO_ACIDS = frozenset("ACDEFGHIKLMNPQRSTVWY")
_MISSING_TEXT = frozenset({"", "NA", "N/A", "NULL", "NONE", "NAN", "<NA>"})


def _is_missing(value: object) -> bool:
    """Return whether a scalar metadata/config value is empty or NA-like."""
    if value is None:
        return True
    if isinstance(value, str) and value.strip().upper() in _MISSING_TEXT:
        return True
    try:
        missing = pd.isna(value)
    except (TypeError, ValueError):
        return False
    if isinstance(missing, bool):
        return missing
    return False


def _require_mapping_fields(record: Mapping[str, Any], required: tuple[str, ...], context: str) -> None:
    missing = [
        field
        for field in required
        if field not in record or _is_missing(record[field])
    ]
    if missing:
        joined = ", ".join(missing)
        raise ValueError(f"{context} is missing required fields: {joined}")


def _normalize_sequence(value: object) -> str:
    if pd.isna(value):
        return ""
    return str(value).strip().upper()


def _canonical_missense(values: pd.Series, sequence: str) -> set[str]:
    mutations: set[str] = set()
    for value in values.dropna().astype(str):
        value = value.strip().upper()
        match = MISSENSE_PATTERN.fullmatch(value)
        if match is None:
            continue
        wild_type, position_text, mutant = match.groups()
        if wild_type not in STANDARD_AMINO_ACIDS or mutant not in STANDARD_AMINO_ACIDS:
            continue
        if wild_type == mutant:
            continue
        position = int(position_text)
        if position < 1 or position > len(sequence) or sequence[position - 1] != wild_type:
            raise ValueError(f"Mutation {value} has a wild-type residue inconsistent with reference")
        mutations.add(value)
    return mutations


def audit_pair_overlap(
    abundance: pd.DataFrame,
    function: pd.DataFrame,
    abundance_sequence: str,
    function_sequence: str,
) -> dict[str, int | bool]:
    """Count exactly shared single substitutions after validating both references."""
    sequence_a = _normalize_sequence(abundance_sequence)
    sequence_f = _normalize_sequence(function_sequence)
    if sequence_a != sequence_f:
        raise ValueError("Paired assay reference sequences differ")
    for name, frame in (("abundance", abundance), ("function", function)):
        if "mutant" not in frame:
            raise ValueError(f"{name} assay is missing required mutant column")

    abundance_mutants = _canonical_missense(abundance["mutant"], sequence_a)
    function_mutants = _canonical_missense(function["mutant"], sequence_a)
    shared = abundance_mutants & function_mutants
    positions = {int(MISSENSE_PATTERN.fullmatch(value).group(2)) for value in shared}
    return {
        "sequence_equal": True,
        "abundance_unique_missense": len(abundance_mutants),
        "function_unique_missense": len(function_mutants),
        "shared_missense": len(shared),
        "shared_positions": len(positions),
        "bad_wild_type": 0,
    }


def build_inventory(reference: pd.DataFrame, candidates: dict) -> pd.DataFrame:
    """Return one provenance-rich row for every configured candidate assay.

    The function validates each selected reference row independently. It does
    not infer a target sequence from another assay with the same protein ID.
    """
    missing_reference_columns = [column for column in REFERENCE_COLUMNS if column not in reference]
    if missing_reference_columns:
        raise ValueError(
            "Reference metadata is missing required columns: "
            + ", ".join(missing_reference_columns)
        )

    source = candidates.get("source")
    assays = candidates.get("assays")
    if not isinstance(source, Mapping) or not isinstance(assays, Mapping):
        raise TypeError("Candidates must define mapping-valued source and assays sections")
    _require_mapping_fields(source, SOURCE_COLUMNS, "Candidates source")

    selected_ids = list(assays)
    if not selected_ids:
        raise ValueError("Candidates must define at least one assay")
    if reference["DMS_id"].duplicated().any():
        duplicates = reference.loc[reference["DMS_id"].duplicated(), "DMS_id"].tolist()
        raise ValueError(f"Reference metadata has duplicate DMS_id values: {duplicates}")

    selected = reference.loc[reference["DMS_id"].isin(selected_ids)].copy()
    missing_ids = sorted(set(selected_ids) - set(selected["DMS_id"]))
    if missing_ids:
        raise ValueError(f"Candidate assay IDs absent from reference metadata: {missing_ids}")

    required_reference_values = (
        ("DMS_filename", "DMS_filename"),
        ("DMS_number_single_mutants", "mutant count"),
        ("UniProt_ID", "UniProt_ID"),
        ("jo", "DOI (jo)"),
        ("taxon", "organism (taxon)"),
        ("source_organism", "source organism"),
        ("selection_assay", "readout (selection_assay)"),
        ("raw_DMS_directionality", "directionality"),
        ("target_seq", "target_seq"),
    )
    missing_values: list[str] = []
    for row in selected.to_dict(orient="records"):
        assay_id = row["DMS_id"]
        for field, label in required_reference_values:
            if _is_missing(row[field]):
                missing_values.append(f"{assay_id}: {label}")
    if missing_values:
        raise ValueError(
            "Selected assay rows missing required reference values: " + "; ".join(missing_values)
        )

    records: list[dict[str, Any]] = []
    for row in selected.to_dict(orient="records"):
        assay_id = row["DMS_id"]
        assay = assays[assay_id]
        if not isinstance(assay, Mapping):
            raise TypeError(f"Candidate assay {assay_id} must be a mapping")
        _require_mapping_fields(assay, CANDIDATE_COLUMNS, f"Candidate assay {assay_id}")

        sequence = _normalize_sequence(row["target_seq"])
        records.append(
            {
                "assay_id": assay_id,
                "reference_filename": row["DMS_filename"],
                "metadata_single_mutant_count": row["DMS_number_single_mutants"],
                "protein_id": row["UniProt_ID"],
                "role": assay["role"],
                "doi": row["jo"],
                "assay_url": assay["assay_url"],
                "metadata_source": source["name"],
                "metadata_url": source["url"],
                "metadata_downloaded_at": source["downloaded_at"],
                "metadata_sha256": source["sha256"],
                "metadata_revision": source["revision"],
                "organism": row["taxon"],
                "source_organism": row["source_organism"],
                "host_system": assay["host_system"],
                "construct": assay["construct"],
                "sequence": sequence,
                "sequence_hash": hashlib.sha256(sequence.encode("ascii")).hexdigest(),
                "sequence_length": len(sequence),
                "tags": json.dumps(assay["tags"], ensure_ascii=True, sort_keys=True),
                "readout": row["selection_assay"],
                "selection_type": row["selection_type"],
                "normalization": assay["normalization"],
                "directionality": row["raw_DMS_directionality"],
                "controls": assay["controls"],
                "qc_rule": assay["qc_rule"],
                "license": assay["license"],
                "threshold_source": assay["threshold_source"],
            }
        )

    return pd.DataFrame.from_records(records)


def validate_pair_config(
    inventory: pd.DataFrame, pairs: Mapping[str, Any]
) -> list[dict[str, str]]:
    """Validate and resolve configured abundance/function assay pairs.

    Returning resolved rows gives callers a safe interface for downstream file
    access and avoids unchecked DataFrame ``.loc`` lookups in reporting scripts.
    """
    if not isinstance(pairs, Mapping) or not pairs:
        raise ValueError("Pair configuration must define at least one pair")
    required_inventory_columns = {"assay_id", "role", "protein_id"}
    missing_columns = sorted(required_inventory_columns - set(inventory.columns))
    if missing_columns:
        raise ValueError(
            "Inventory is missing required pair-validation columns: " + ", ".join(missing_columns)
        )
    if inventory["assay_id"].duplicated().any():
        duplicates = inventory.loc[inventory["assay_id"].duplicated(), "assay_id"].tolist()
        raise ValueError(f"Inventory has duplicate assay IDs: {duplicates}")

    by_id = {
        row["assay_id"]: row
        for row in inventory.to_dict(orient="records")
        if not _is_missing(row["assay_id"])
    }
    if len(by_id) != len(inventory):
        raise ValueError("Inventory contains missing assay IDs")

    resolved: list[dict[str, str]] = []
    definitions: set[tuple[str, str]] = set()
    used_assays: dict[str, str] = {}
    for pair_name, pair in pairs.items():
        if _is_missing(pair_name):
            raise ValueError("Pair definition has a missing name")
        if not isinstance(pair, Mapping):
            raise TypeError(f"Pair {pair_name} must be a mapping")
        unknown_keys = sorted(set(pair) - {"abundance", "function"})
        if unknown_keys:
            raise ValueError(f"Pair {pair_name} has invalid fields: {unknown_keys}")
        missing_roles = [role for role in ("abundance", "function") if _is_missing(pair.get(role))]
        if missing_roles:
            raise ValueError(
                f"Pair {pair_name} is missing {missing_roles[0]} assay ID"
            )
        abundance_id = pair["abundance"]
        function_id = pair["function"]
        if not isinstance(abundance_id, str) or not isinstance(function_id, str):
            raise TypeError(f"Pair {pair_name} assay IDs must be strings")
        definition = (abundance_id, function_id)
        if definition in definitions:
            raise ValueError(f"Pair {pair_name} has duplicate pair definition {definition}")
        definitions.add(definition)
        if abundance_id == function_id:
            raise ValueError(f"Pair {pair_name} uses the same assay for abundance and function")
        for role, assay_id in (("abundance", abundance_id), ("function", function_id)):
            if assay_id not in by_id:
                raise ValueError(f"Pair {pair_name} references unknown {role} assay ID: {assay_id}")
            previous_pair = used_assays.get(assay_id)
            if previous_pair is not None:
                raise ValueError(
                    f"Pair {pair_name} reuses assay ID {assay_id} already used by {previous_pair}"
                )
            used_assays[assay_id] = str(pair_name)
            row = by_id[assay_id]
            if row["role"] != role:
                raise ValueError(
                    f"Pair {pair_name} {role} assay {assay_id} has role {row['role']!r}; "
                    f"expected {role!r}"
                )

        abundance_row = by_id[abundance_id]
        function_row = by_id[function_id]
        if _is_missing(abundance_row["protein_id"]) or _is_missing(function_row["protein_id"]):
            raise ValueError(f"Pair {pair_name} has missing protein ID")
        if abundance_row["protein_id"] != function_row["protein_id"]:
            raise ValueError(
                f"Pair {pair_name} has protein mismatch: "
                f"{abundance_row['protein_id']} != {function_row['protein_id']}"
            )
        resolved.append(
            {
                "pair_name": str(pair_name),
                "abundance_assay": abundance_id,
                "function_assay": function_id,
                "protein_id": str(abundance_row["protein_id"]),
            }
        )
    return resolved
