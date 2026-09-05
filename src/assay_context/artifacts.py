"""Public-tree policy for row-level, source-derived analysis artifacts."""

from __future__ import annotations

from pathlib import Path

EXTERNAL_ARTIFACT_DIR = Path("D:/protein_assay_work/release_artifacts")
ROW_LEVEL_ARTIFACTS = frozenset(
    {
        "assay_atlas.csv",
        "task4_oof_predictions.csv",
        "task4_lopo_predictions.csv",
        "task5_oof_predictions.csv",
        "task5_lopo_predictions.csv",
        "task5_complete_case_oof_predictions.csv",
        "task5_structure_mapping.csv",
    }
)


def external_artifact_path(filename: str) -> Path:
    """Return the D-drive path for one approved row-level output filename."""
    path = Path(filename)
    if path.name != filename or filename not in ROW_LEVEL_ARTIFACTS:
        raise ValueError(f"unsupported row-level artifact filename: {filename!r}")
    return EXTERNAL_ARTIFACT_DIR / filename
