from __future__ import annotations

import hashlib
import importlib.util
from pathlib import Path

import pandas as pd
import pytest

from assay_context.inventory import audit_pair_overlap, build_inventory, validate_pair_config

_SCRIPT_SPEC = importlib.util.spec_from_file_location("inventory_script", Path("scripts/01_inventory.py"))
assert _SCRIPT_SPEC and _SCRIPT_SPEC.loader
_01_inventory = importlib.util.module_from_spec(_SCRIPT_SPEC)
_SCRIPT_SPEC.loader.exec_module(_01_inventory)


@pytest.fixture
def candidates() -> dict:
    return {
        "source": {
            "name": "ProteinGym DMS substitutions metadata",
            "url": "https://example.test/metadata.csv",
            "downloaded_at": "2026-09-05",
            "sha256": "a" * 64,
            "revision": "test-revision",
        },
        "assays": {
            "ABUNDANCE_A": {
                "role": "abundance",
                "assay_url": "https://doi.org/10.1000/example",
                "construct": "full-length protein fused to GFP",
                "host_system": "HEK293T cells",
                "tags": ["GFP"],
                "normalization": "author-reported normalization",
                "controls": "wild type and nonsense controls",
                "qc_rule": "author-reported quality filters",
                "threshold_source": "author methods",
                "license": "source repository terms",
            },
            "FUNCTION_A": {
                "role": "function",
                "assay_url": "https://doi.org/10.1000/example",
                "construct": "full-length protein",
                "host_system": "yeast complementation",
                "tags": [],
                "normalization": "author-reported normalization",
                "controls": "wild type and empty-vector controls",
                "qc_rule": "author-reported quality filters",
                "threshold_source": "author methods",
                "license": "source repository terms",
            },
        },
    }


@pytest.fixture
def paired_reference() -> pd.DataFrame:
    return pd.DataFrame(
        [
            {
                "DMS_id": "ABUNDANCE_A",
                "DMS_filename": "ABUNDANCE_A.csv",
                "DMS_number_single_mutants": 100,
                "UniProt_ID": "TEST_HUMAN",
                "target_seq": "ACDE",
                "jo": "10.1000/example",
                "taxon": "Human",
                "source_organism": "Homo sapiens",
                "selection_assay": "protein abundance",
                "selection_type": "Fluorescence",
                "raw_DMS_directionality": "1",
                "DMS_binarization_method": "manual",
                "DMS_binarization_cutoff": "0.5",
            },
            {
                "DMS_id": "FUNCTION_A",
                "DMS_filename": "FUNCTION_A.csv",
                "DMS_number_single_mutants": 101,
                "UniProt_ID": "TEST_HUMAN",
                "target_seq": "ACDE",
                "jo": "10.1000/example",
                "taxon": "Human",
                "source_organism": "Homo sapiens",
                "selection_assay": "growth complementation",
                "selection_type": "Activity",
                "raw_DMS_directionality": "1",
                "DMS_binarization_method": "manual",
                "DMS_binarization_cutoff": "0.4",
            },
        ]
    )


def test_build_inventory_records_one_row_per_candidate_assay(
    paired_reference: pd.DataFrame, candidates: dict
) -> None:
    """A missing candidate or an incorrect sequence hash must not silently enter the audit."""
    inventory = build_inventory(paired_reference, candidates)

    assert inventory["assay_id"].tolist() == ["ABUNDANCE_A", "FUNCTION_A"]
    assert inventory["protein_id"].tolist() == ["TEST_HUMAN", "TEST_HUMAN"]
    assert inventory["sequence_hash"].tolist() == [hashlib.sha256(b"ACDE").hexdigest()] * 2
    assert inventory["role"].tolist() == ["abundance", "function"]
    assert inventory["reference_filename"].tolist() == ["ABUNDANCE_A.csv", "FUNCTION_A.csv"]
    assert inventory["metadata_single_mutant_count"].tolist() == [100, 101]
    assert inventory["metadata_url"].tolist() == ["https://example.test/metadata.csv"] * 2
    assert inventory["metadata_sha256"].tolist() == ["a" * 64] * 2


def test_build_inventory_rejects_assay_row_without_target_sequence(
    paired_reference: pd.DataFrame, candidates: dict
) -> None:
    """Removing target_seq validation would make an unsafe name-only assay pairing possible."""
    paired_reference.loc[1, "target_seq"] = ""

    with pytest.raises(ValueError, match="FUNCTION_A.*target_seq"):
        build_inventory(paired_reference, candidates)


@pytest.mark.parametrize(
    ("field", "message"),
    [
        ("DMS_filename", "DMS_filename"),
        ("DMS_number_single_mutants", "mutant count"),
        ("UniProt_ID", "UniProt_ID"),
        ("jo", "DOI"),
        ("taxon", "organism"),
        ("selection_assay", "readout"),
        ("raw_DMS_directionality", "directionality"),
    ],
)
def test_build_inventory_rejects_missing_required_reference_values(
    paired_reference: pd.DataFrame, candidates: dict, field: str, message: str
) -> None:
    paired_reference.loc[0, field] = pd.NA

    with pytest.raises(ValueError, match=f"ABUNDANCE_A.*{message}"):
        build_inventory(paired_reference, candidates)


def test_validate_pair_config_returns_resolved_rows(
    paired_reference: pd.DataFrame, candidates: dict
) -> None:
    inventory = build_inventory(paired_reference, candidates)

    result = validate_pair_config(
        inventory,
        {"TEST_HUMAN": {"abundance": "ABUNDANCE_A", "function": "FUNCTION_A"}},
    )

    assert result == [
        {
            "pair_name": "TEST_HUMAN",
            "abundance_assay": "ABUNDANCE_A",
            "function_assay": "FUNCTION_A",
            "protein_id": "TEST_HUMAN",
        }
    ]


def test_overlap_output_row_uses_configured_pair_name() -> None:
    pair = {
        "pair_name": "CYP2C9",
        "protein_id": "CP2C9_HUMAN",
        "abundance_assay": "ABUNDANCE_A",
        "function_assay": "FUNCTION_A",
    }

    row = _01_inventory.overlap_output_row(
        pair,
        {"shared_missense": 400},
        mapping_gate="PASS",
        mapping_gate_reason="meets thresholds",
        abundance_file_sha256="a" * 64,
        function_file_sha256="b" * 64,
    )

    assert row["protein"] == "CYP2C9"
    assert row["abundance_assay"] == "ABUNDANCE_A"
    assert row["function_assay"] == "FUNCTION_A"


def test_unique_reasons_preserves_first_seen_order() -> None:
    assert _01_inventory.unique_reasons(["same", "same", "different", "same"]) == [
        "same",
        "different",
    ]


def test_validate_pair_config_rejects_missing_assay_id(
    paired_reference: pd.DataFrame, candidates: dict
) -> None:
    inventory = build_inventory(paired_reference, candidates)

    with pytest.raises(ValueError, match="missing abundance assay ID"):
        validate_pair_config(inventory, {"TEST_HUMAN": {"function": "FUNCTION_A"}})


def test_validate_pair_config_rejects_role_mismatch(
    paired_reference: pd.DataFrame, candidates: dict
) -> None:
    inventory = build_inventory(paired_reference, candidates)

    with pytest.raises(ValueError, match="abundance assay .* role"):
        validate_pair_config(
            inventory,
            {"TEST_HUMAN": {"abundance": "FUNCTION_A", "function": "ABUNDANCE_A"}},
        )


def test_validate_pair_config_rejects_protein_mismatch(
    paired_reference: pd.DataFrame, candidates: dict
) -> None:
    paired_reference.loc[1, "UniProt_ID"] = "OTHER_HUMAN"
    inventory = build_inventory(paired_reference, candidates)

    with pytest.raises(ValueError, match="protein mismatch"):
        validate_pair_config(
            inventory,
            {"TEST_HUMAN": {"abundance": "ABUNDANCE_A", "function": "FUNCTION_A"}},
        )


def test_validate_pair_config_rejects_duplicate_definition(
    paired_reference: pd.DataFrame, candidates: dict
) -> None:
    inventory = build_inventory(paired_reference, candidates)

    with pytest.raises(ValueError, match="duplicate pair definition"):
        validate_pair_config(
            inventory,
            {
                "FIRST": {"abundance": "ABUNDANCE_A", "function": "FUNCTION_A"},
                "SECOND": {"abundance": "ABUNDANCE_A", "function": "FUNCTION_A"},
            },
        )


def test_audit_pair_overlap_counts_shared_missense_and_validates_wild_type() -> None:
    sequence = "ACDE"
    abundance = pd.DataFrame(
        {"mutant": ["A1V", "C2A", "D3*", "A1V:C2A", "A1A", "B1V", "C2U"]}
    )
    function = pd.DataFrame({"mutant": ["A1V", "C2A", "E4K", "A1A", "B1V"]})

    result = audit_pair_overlap(abundance, function, sequence, sequence)

    assert result == {
        "sequence_equal": True,
        "abundance_unique_missense": 2,
        "function_unique_missense": 3,
        "shared_missense": 2,
        "shared_positions": 2,
        "bad_wild_type": 0,
    }


def test_audit_pair_overlap_rejects_reference_mismatch() -> None:
    abundance = pd.DataFrame({"mutant": ["A1V"]})
    function = pd.DataFrame({"mutant": ["A1V"]})

    with pytest.raises(ValueError, match="reference sequences differ"):
        audit_pair_overlap(abundance, function, "ACDE", "ACDF")


def test_audit_pair_overlap_rejects_wrong_wild_type() -> None:
    abundance = pd.DataFrame({"mutant": ["G1V"]})
    function = pd.DataFrame({"mutant": ["G1V"]})

    with pytest.raises(ValueError, match="wild-type residue"):
        audit_pair_overlap(abundance, function, "ACDE", "ACDE")
