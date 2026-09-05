from __future__ import annotations

import pandas as pd

from assay_context.atlas import summarize_phenotypes


def test_summarize_phenotypes_has_explicit_denominator_and_expected_categories() -> None:
    pairs = pd.DataFrame(
        {
            "protein": ["TEST"] * 5,
            "variant": ["A1V", "A1L", "A1I", "A1F", "A1Y"],
            "position": [1, 1, 2, 2, 3],
            "abundance_label": ["preserved", "low", "preserved", "low", "uncertain"],
            "function_label": ["impaired", "impaired", "preserved", "preserved", "uncertain"],
            "primary_endpoint": [
                "function-impaired",
                "not_eligible",
                "function-preserved",
                "not_eligible",
                "not_eligible",
            ],
        }
    )

    result = summarize_phenotypes(pairs)

    assert set(result["phenotype"]) == {
        "abundance-preserved_function-impaired",
        "abundance-preserved_function-preserved",
        "abundance-low_function-impaired",
        "abundance-low_function-preserved",
        "uncertain-or-other",
    }
    assert result["denominator"].eq(5).all()
    assert result["count"].sum() == 5
    assert result["fraction"].sum() == 1.0


def test_summarize_phenotypes_keeps_proteins_and_zero_count_categories() -> None:
    pairs = pd.DataFrame(
        {
            "protein": ["A", "B"],
            "variant": ["A1V", "A1V"],
            "position": [1, 1],
            "abundance_label": ["preserved", "preserved"],
            "function_label": ["impaired", "impaired"],
            "primary_endpoint": ["function-impaired"] * 2,
        }
    )

    result = summarize_phenotypes(pairs)

    assert set(result["protein"]) == {"A", "B"}
    assert result.groupby("protein").size().eq(5).all()
    assert result.loc[result["phenotype"].eq("abundance-preserved_function-impaired"), "count"].eq(1).all()
    assert result.loc[result["phenotype"].eq("abundance-low_function-preserved"), "count"].eq(0).all()
