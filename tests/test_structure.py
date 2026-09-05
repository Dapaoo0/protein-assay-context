from __future__ import annotations

import importlib.util
from pathlib import Path

import numpy as np
import pytest
import yaml

from assay_context.structure import map_structure

_SCRIPT_SPEC = importlib.util.spec_from_file_location(
    "structure_sensitivity_script", Path("scripts/05_structure_sensitivity.py")
)
assert _SCRIPT_SPEC and _SCRIPT_SPEC.loader
_05_structure_sensitivity = importlib.util.module_from_spec(_SCRIPT_SPEC)
_SCRIPT_SPEC.loader.exec_module(_05_structure_sensitivity)


def _pdb(path: Path) -> Path:
    # Residue 2 is intentionally absent from the observed chain.  Residue 3
    # has an insertion code so alignment is exercised rather than residue
    # numbering being treated as a sufficient mapping.
    lines = [
        "HEADER    TEST STRUCTURE 1TST",
        "HELIX    1   1 ALA A   1  ALA A   1  1",
        "ATOM      1  N   ALA A   1      0.000   0.000   0.000  1.00 11.00           N",
        "ATOM      2  CA  ALA A   1      1.000   0.000   0.000  1.00 12.00           C",
        "ATOM      3  C   ALA A   1      1.500   1.000   0.000  1.00 13.00           C",
        "ATOM      4  O   ALA A   1      1.000   2.000   0.000  1.00 14.00           O",
        "ATOM      5  N   GLY A   3      4.000   0.000   0.000  1.00 21.00           N",
        "ATOM      6  CA  GLY A   3      5.000   0.000   0.000  1.00 22.00           C",
        "ATOM      7  C   GLY A   3      5.500   1.000   0.000  1.00 23.00           C",
        "ATOM      8  O   GLY A   3      5.000   2.000   0.000  1.00 24.00           O",
        "TER",
        "END",
    ]
    path.write_text("\n".join(lines) + "\n", encoding="ascii")
    return path


def test_map_structure_emits_every_reference_position_and_missing_is_nan(tmp_path: Path) -> None:
    result = map_structure("AGG", str(_pdb(tmp_path / "1tst.pdb")), "A")
    assert result["position"].tolist() == [1, 2, 3]
    missing = result.loc[result["position"].eq(2)].iloc[0]
    assert missing["mapping_status"] == "missing"
    assert np.isnan(missing["accessibility"])
    assert np.isnan(missing["contacts"])
    assert np.isnan(missing["site_distance"])
    assert result.loc[result["position"].eq(1), "mapping_status"].item() == "mapped"


def test_experimental_b_factor_is_not_reported_as_plddt(tmp_path: Path) -> None:
    result = map_structure("AGG", str(_pdb(tmp_path / "1tst.pdb")), "A")
    assert result["confidence_type"].dropna().unique().tolist() == ["B-factor"]
    assert result["confidence_value"].notna().any()
    assert result["plddt"].isna().all()


def test_wrong_chain_and_low_identity_fail_closed(tmp_path: Path) -> None:
    path = _pdb(tmp_path / "1tst.pdb")
    with pytest.raises(ValueError, match="chain"):
        map_structure("AGG", str(path), "B")
    with pytest.raises(ValueError, match="identity|coverage"):
        map_structure("WWW", str(path), "A")


def test_functional_site_distance_uses_predeclared_positions(tmp_path: Path) -> None:
    result = map_structure(
        "AGG", str(_pdb(tmp_path / "1tst.pdb")), "A", functional_positions=[1],
        functional_site_source="test annotation",
    )
    assert result.loc[result["position"].eq(1), "site_distance"].notna().item()
    assert result.attrs["functional_site_source"] == "test annotation"


def test_structure_stage_restores_release_freeze_without_dropping_prior_sections(
    tmp_path: Path,
) -> None:
    config = tmp_path / "analysis_freeze.yaml"
    config.write_text(
        yaml.safe_dump({"model_evaluation": {"gate_b": "PASS"}}, sort_keys=False),
        encoding="utf-8",
    )

    _05_structure_sensitivity.update_analysis_freeze(config)

    freeze = yaml.safe_load(config.read_text(encoding="utf-8"))
    assert freeze["model_evaluation"] == {"gate_b": "PASS"}
    structure = freeze["structure_evaluation"]
    assert structure["primary_proteins"] == ["CYP2C9", "PTEN"]
    assert structure["structures"]["CYP2C9"]["pdb_id"] == "1R9O"
    assert structure["structures"]["PTEN"]["functional_positions"] == [124]
    assert structure["case_rule"].startswith(
        "deterministically selected post-analysis using a recorded rule"
    )
