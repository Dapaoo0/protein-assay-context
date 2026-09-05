from __future__ import annotations

import importlib.util
from pathlib import Path
from zipfile import ZipFile

import pandas as pd
import pytest

from assay_context.scores import audit_cross_assay_scores, join_scores

_SCRIPT_SPEC = importlib.util.spec_from_file_location("build_atlas", Path("scripts/03_build_atlas.py"))
assert _SCRIPT_SPEC and _SCRIPT_SPEC.loader
_03_build_atlas = importlib.util.module_from_spec(_SCRIPT_SPEC)
_SCRIPT_SPEC.loader.exec_module(_03_build_atlas)


def _pairs() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "protein": ["TEST", "TEST", "TEST"],
            "protein_id": ["TEST_HUMAN"] * 3,
            "sequence_hash": ["a" * 64] * 3,
            "variant": ["A1V", "A1L", "A1I"],
            "wt": ["A"] * 3,
            "position": [1] * 3,
            "mutant": ["V", "L", "I"],
        }
    )


def test_missing_model_scores_preserve_every_paired_variant() -> None:
    scores = pd.DataFrame(
        {
            "protein": ["TEST", "TEST"],
            "protein_id": ["TEST_HUMAN"] * 2,
            "sequence_hash": ["a" * 64] * 2,
            "wt": ["A"] * 2,
            "position": [1] * 2,
            "mutant": ["V", "I"],
            "variant": ["A1V", "A1I"],
            "model": ["ESM1v_single"] * 2,
            "score": [0.8, 0.2],
        }
    )

    result = join_scores(_pairs(), scores)

    assert result["variant"].tolist() == ["A1V", "A1L", "A1I"]
    assert result.loc[0, "ESM1v_single"] == 0.8
    assert pd.isna(result.loc[1, "ESM1v_single"])
    assert result.loc[2, "ESM1v_single"] == 0.2
    assert result.attrs["coverage"]["ESM1v_single"]["matched"] == 2
    assert result.attrs["coverage"]["ESM1v_single"]["total"] == 3


def test_duplicate_canonical_score_key_is_rejected_even_when_values_match() -> None:
    scores = pd.DataFrame(
        {
            "protein": ["TEST", "TEST"],
            "protein_id": ["TEST_HUMAN"] * 2,
            "sequence_hash": ["a" * 64] * 2,
            "wt": ["A"] * 2,
            "position": [1] * 2,
            "mutant": ["V"] * 2,
            "variant": ["A1V", "A1V"],
            "assay_id": ["TEST_abundance", "TEST_function"],
            "model": ["ESM1v_single"] * 2,
            "score": [0.8, 0.8],
        }
    )

    with pytest.raises(ValueError, match="duplicate canonical score keys"):
        join_scores(_pairs(), scores)


def test_same_model_variant_mismatch_is_recorded_and_blocks_merge() -> None:
    scores = pd.DataFrame(
        {
            "protein": ["TEST", "TEST"],
            "protein_id": ["TEST_HUMAN"] * 2,
            "sequence_hash": ["a" * 64] * 2,
            "wt": ["A"] * 2,
            "position": [1] * 2,
            "mutant": ["V"] * 2,
            "variant": ["A1V", "A1V"],
            "assay_id": ["TEST_abundance", "TEST_function"],
            "model": ["ESM1v_single"] * 2,
            "score": [0.8, 0.7],
        }
    )

    with pytest.raises(ValueError, match="duplicate canonical score keys"):
        join_scores(_pairs(), scores)


def test_cross_assay_audit_reports_nonidentical_values_without_coalescing() -> None:
    abundance = _pairs().iloc[[0]].assign(
        assay_id="TEST_abundance", model="ESM1v_single", score=0.8
    )
    function = _pairs().iloc[[0]].assign(
        assay_id="TEST_function", model="ESM1v_single", score=0.80001
    )

    summary, detail = audit_cross_assay_scores(abundance, function)

    assert summary == {
        "groups_checked": 1,
        "exact_match_count": 0,
        "nonidentical_count": 1,
        "abundance_only_count": 0,
        "function_only_count": 0,
    }
    assert detail.loc[0, "assay_id_1"] == "TEST_abundance"
    assert detail.loc[0, "assay_id_2"] == "TEST_function"
    assert detail.loc[0, "absolute_difference"] == pytest.approx(0.00001)


def test_variant_token_must_match_canonical_components() -> None:
    pairs = _pairs().copy()
    pairs.loc[0, "variant"] = "A2V"

    with pytest.raises(ValueError, match="inconsistent variant"):
        join_scores(pairs, pd.DataFrame(columns=["protein", "protein_id", "sequence_hash", "variant", "wt", "position", "mutant", "model", "score"]))


def test_sequence_hash_is_part_of_score_join_key() -> None:
    scores = pd.DataFrame(
        {
            "protein": ["TEST"],
            "protein_id": ["TEST_HUMAN"],
            "sequence_hash": ["b" * 64],
            "variant": ["A1V"],
            "wt": ["A"],
            "position": [1],
            "mutant": ["V"],
            "model": ["ESM1v_single"],
            "score": [0.8],
        }
    )
    with pytest.raises(ValueError, match="sequence_hash mismatch"):
        join_scores(_pairs(), scores)


@pytest.mark.parametrize("column", ["construct", "offset"])
def test_construct_and_offset_are_validated_when_present(column: str) -> None:
    pairs = _pairs()
    pairs[column] = "reference" if column == "construct" else 0
    scores = pd.DataFrame(
        {
            "protein": ["TEST"],
            "protein_id": ["TEST_HUMAN"],
            "sequence_hash": ["a" * 64],
            "variant": ["A1V"],
            "wt": ["A"],
            "position": [1],
            "mutant": ["V"],
            "model": ["ESM1v_single"],
            "score": [0.8],
            column: ["different" if column == "construct" else 1],
        }
    )

    with pytest.raises(ValueError, match=f"{column} mismatch"):
        join_scores(pairs, scores)


def test_changed_extracted_score_file_fails_against_frozen_hash(tmp_path) -> None:
    archive = tmp_path / "scores.zip"
    score_dir = tmp_path / "scores"
    filename = next(iter(_03_build_atlas.ASSAYS))
    with ZipFile(archive, "w") as zipped:
        zipped.writestr(f"DMS_ProteinGym_substitutions/{filename}", "immutable")
    score_dir.mkdir()
    (score_dir / filename).write_text("changed", encoding="utf-8")

    with pytest.raises(ValueError, match="Frozen extracted score SHA256 mismatch"):
        _03_build_atlas.ensure_extracted_scores(score_dir, archive, {filename: "0" * 64})


def test_identical_metadata_generation_preserves_download_event_timestamp() -> None:
    prior = {"score_archive": {"sha256": "abc", "downloaded_at": "2026-09-05T10:30:13+00:00"}}

    first = _03_build_atlas.frozen_download_timestamp(prior, "abc", "2026-09-06T00:00:00+00:00")
    second = _03_build_atlas.frozen_download_timestamp(prior, "abc", "2026-09-07T00:00:00+00:00")

    assert first == second == "2026-09-05T10:30:13+00:00"


def test_descriptive_flow_label_names_missing_primary_benchmark_rows() -> None:
    pairs = _pairs().iloc[[0]].assign(protein="GCK", primary_eligible=False, endpoint_labeled=False)
    joined = pairs.assign(ESM1v_single=pd.NA, ESM2_650M=pd.NA)
    attrition = pd.DataFrame(
        {
            "protein": ["GCK"],
            "abundance_input": [10],
            "function_input": [11],
            "shared_proteingym_pairs": [1],
            "author_annotated_pairs": [1],
        }
    )

    flow = _03_build_atlas.cohort_flow_table(pairs, joined, attrition)

    assert flow["funnel"].unique().tolist() == ["descriptive-only (no primary benchmark rows)"]


def test_cohort_flow_separates_endpoint_labeled_from_primary_eligible() -> None:
    pairs = _pairs().iloc[[0]].assign(
        protein="GCK", primary_eligible=False, endpoint_labeled=True,
        primary_endpoint="function-impaired",
    )
    joined = pairs.assign(ESM1v_single=pd.NA, ESM2_650M=pd.NA)
    attrition = pd.DataFrame(
        {
            "protein": ["GCK"],
            "abundance_input": [10],
            "function_input": [11],
            "shared_proteingym_pairs": [1],
            "author_annotated_pairs": [1],
        }
    )

    flow = _03_build_atlas.cohort_flow_table(pairs, joined, attrition)

    counts = flow.set_index("stage")["count"]
    assert counts["endpoint labeled"] == 1
    assert counts["primary benchmark eligible"] == 0
