from __future__ import annotations

import pandas as pd
import pytest

from assay_context.variants import pair_assays, parse_variant


def test_parse_variant_validates_reference_sequence() -> None:
    assert parse_variant("A1V", "ACDE") == ("A", 1, "V")


@pytest.mark.parametrize("text", ["G1V", "A0V", "A5V", "A1V:C2A", "A1*", "A1="])
def test_parse_variant_rejects_wrong_or_non_missense_variant(text: str) -> None:
    with pytest.raises(ValueError):
        parse_variant(text, "ACDE")


def test_parse_variant_rejects_non_string_inputs() -> None:
    with pytest.raises(TypeError):
        parse_variant(None, "ACDE")  # type: ignore[arg-type]


def test_pair_assays_is_one_to_one_and_rejects_duplicate_canonical_keys() -> None:
    abundance = pd.DataFrame({"mutant": ["A1V", "C2A"], "DMS_score": [0.8, 0.7]})
    function = pd.DataFrame({"variant": ["A1V", "C2A"], "DMS_score": [0.2, 0.9]})
    paired = pair_assays(abundance, function)
    assert paired[["wt", "position", "mutant"]].to_dict("records") == [
        {"wt": "A", "position": 1, "mutant": "V"},
        {"wt": "C", "position": 2, "mutant": "A"},
    ]
    assert len(paired) == 2

    duplicate = pd.concat([abundance, abundance.iloc[[0]]], ignore_index=True)
    with pytest.raises(ValueError, match="duplicate canonical"):
        pair_assays(duplicate, function)


def test_pair_assays_exposes_attrition_audit() -> None:
    paired = pair_assays(
        pd.DataFrame({"mutant": ["A1V", "C2A"]}),
        pd.DataFrame({"mutant": ["A1V", "D3A"]}),
    )
    assert len(paired) == 1
    assert paired.attrs["attrition"] == {
        "abundance_input": 2,
        "function_input": 2,
        "abundance_valid_single": 2,
        "function_valid_single": 2,
        "shared_pairs": 1,
        "abundance_unmatched": 1,
        "function_unmatched": 1,
    }
