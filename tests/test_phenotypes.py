from __future__ import annotations

import pandas as pd
import pytest

from assay_context.phenotypes import (
    adapt_cyp2c9,
    adapt_gck,
    adapt_pten_abundance,
    adapt_pten_function,
    adapt_vkor,
    gate_a_status,
    label_interval,
    primary_binary_eligibility,
)


def test_label_interval_strict_boundaries() -> None:
    assert label_interval(0.1, 0.4, 0.5) == "low"
    assert label_interval(0.6, 0.8, 0.5) == "preserved"
    assert label_interval(0.4, 0.6, 0.5) == "uncertain"
    assert label_interval(0.5, 0.7, 0.5) == "uncertain"


def test_cyp_adapter_uses_author_classes_and_point_estimate_flag() -> None:
    source = pd.DataFrame(
        {
            "variant": ["A1V", "A1L", "A1I", "A1F"],
            "class": ["missense"] * 4,
            "activity_class": ["decreased", "wt-like", "possibly_decreased", "increased"],
            "abundance_class": ["wt-like", "decreased", "wt-like", "wt-like"],
        }
    )
    result = adapt_cyp2c9(source)
    assert result[["variant", "abundance_label", "function_label"]].to_dict("records") == [
        {"variant": "A1V", "abundance_label": "preserved", "function_label": "impaired"},
        {"variant": "A1L", "abundance_label": "low", "function_label": "preserved"},
        {"variant": "A1I", "abundance_label": "preserved", "function_label": "uncertain"},
        {"variant": "A1F", "abundance_label": "preserved", "function_label": "other"},
    ]
    assert result["uncertainty_kind"].eq("point-estimate author category").all()


def test_cyp_nonsense_like_function_class_is_impaired() -> None:
    result = adapt_cyp2c9(
        pd.DataFrame(
            {
                "variant": ["A1V"],
                "class": ["missense"],
                "activity_class": ["nonsense-like"],
                "abundance_class": ["wt-like"],
            }
        )
    )
    assert result.loc[0, "function_label"] == "impaired"


def test_cyp_nonsense_like_abundance_class_is_low() -> None:
    result = adapt_cyp2c9(
        pd.DataFrame(
            {
                "variant": ["A1V"],
                "class": ["missense"],
                "activity_class": ["wt-like"],
                "abundance_class": ["nonsense-like"],
            }
        )
    )
    assert result.loc[0, "abundance_label"] == "low"


def test_pten_function_requires_high_confidence_and_uses_paper_cutoffs() -> None:
    source = pd.DataFrame(
        {
            "Variant (one letter)": ["A1V", "A1L", "A1I", "A1F", "A1Y"],
            "Type": ["missense"] * 5,
            "Cum_score": [-1.11, -1.1101, 0.89, 0.8901, -2.0],
            "High_conf": [True, True, True, True, False],
        }
    )
    result = adapt_pten_function(source)
    assert result["function_label"].tolist() == ["impaired", "impaired", "preserved", "other", "uncertain"]


def test_gck_ci_crossings_are_uncertain_and_threshold_equality_is_not_certain() -> None:
    source = pd.DataFrame(
        {
            "variant": ["A1V", "A1L", "A1I"],
            "abundance_score": [0.6, 0.7, 0.5],
            "abundance_score_se": [0.0, 0.0, 0.0],
            "df": [20, 20, 20],
        }
    )
    result = adapt_gck(source, role="abundance")
    assert result["abundance_label"].tolist() == ["uncertain", "preserved", "low"]
    assert result["uncertainty_kind"].str.startswith("95% CI").all()


def test_gck_unexpected_df_fails_explicitly() -> None:
    source = pd.DataFrame(
        {"variant": ["A1V"], "abundance_score": [0.7], "abundance_score_se": [0.1], "df": [25]}
    )
    with pytest.raises(ValueError, match="unsupported GCK degrees of freedom"):
        adapt_gck(source, role="abundance")


def test_pten_and_gck_adapters_preserve_metadata() -> None:
    pten = adapt_pten_abundance(
        pd.DataFrame({"variant": ["A1V"], "abundance_class": ["wt-like"]})
    )
    gck = adapt_gck(
        pd.DataFrame({"variant": ["A1V"], "abundance_score": [0.7], "abundance_score_se": [0.0], "df": [20]}),
        role="abundance",
    )
    assert pten.loc[0, "uncertainty_kind"] and pten.loc[0, "source_adapter"]
    assert gck.loc[0, "uncertainty_kind"] and gck.loc[0, "source_adapter"]


def test_vkor_adapter_keeps_high_and_possible_out_of_preserved() -> None:
    source = pd.DataFrame(
        {
            "variant": ["A1V", "A1L", "A1I", "A1F"],
            "class": ["missense"] * 4,
            "abundance_class": ["wt-like", "low", "possibly_wt-like", "high"],
            "activity_class": ["wt-like", "low", "possibly_wt-like", "high"],
        }
    )
    result = adapt_vkor(source)
    assert result["abundance_label"].tolist() == ["preserved", "low", "uncertain", "other"]
    assert result["function_label"].tolist() == ["preserved", "impaired", "uncertain", "other"]


def test_adapters_reject_duplicate_variant_rows() -> None:
    with pytest.raises(ValueError, match="duplicate"):
        adapt_pten_abundance(
            pd.DataFrame({"variant": ["A1V", "A1V"], "abundance_class": ["wt-like", "low"]})
        )


def test_label_interval_handles_pd_na() -> None:
    assert label_interval(pd.NA, 0.4, 0.5) == "uncertain"


def test_gate_a_reduces_scope_when_total_pairs_or_positions_are_insufficient() -> None:
    frame = pd.DataFrame(
        {
            "position": [1] * 50 + [2] * 50,
            "primary_endpoint": ["function-impaired"] * 50 + ["function-preserved"] * 50,
        }
    )
    assert gate_a_status(frame) == "REDUCE_SCOPE"


def test_primary_binary_eligibility_requires_gate_and_binary_row_label() -> None:
    frame = pd.DataFrame(
        {"primary_endpoint": ["function-impaired", "function-preserved", "not_eligible"]}
    )
    assert primary_binary_eligibility(frame, "PASS").tolist() == [True, True, False]
    assert primary_binary_eligibility(frame, "REDUCE_SCOPE").tolist() == [False, False, False]


def test_primary_benchmark_protein_gate_policy() -> None:
    frame = pd.DataFrame(
        {"primary_endpoint": ["function-impaired", "function-preserved"]}
    )
    assert primary_binary_eligibility(frame, "PASS").all()  # CYP2C9/PTEN policy
    assert not primary_binary_eligibility(frame, "REDUCE_SCOPE").any()  # GCK/VKORC1 policy
