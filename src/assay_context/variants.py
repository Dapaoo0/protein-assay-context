"""Strict variant parsing and one-to-one assay pairing."""

from __future__ import annotations

import re

import pandas as pd

_VARIANT_RE = re.compile(r"^([A-Z])([1-9][0-9]*)([A-Z])$")
_AMINO_ACIDS = frozenset("ACDEFGHIKLMNPQRSTVWY")


def parse_variant(text: str, sequence: str) -> tuple[str, int, str]:
    """Parse and validate one one-letter missense substitution.

    Positions are one-based, as in ProteinGym and the author tables.  The
    sequence check is deliberately strict so that an apparently valid token
    cannot silently be paired against a different construct or offset.
    """
    if not isinstance(text, str) or not isinstance(sequence, str):
        raise TypeError("variant and sequence must be strings")
    token = text.strip().upper()
    reference = sequence.strip().upper()
    match = _VARIANT_RE.fullmatch(token)
    if match is None:
        raise ValueError(f"Not a single amino-acid missense variant: {text!r}")
    wild_type, position_text, mutant = match.groups()
    if wild_type not in _AMINO_ACIDS or mutant not in _AMINO_ACIDS:
        raise ValueError(f"Variant contains a non-standard amino acid: {text!r}")
    position = int(position_text)
    if position > len(reference):
        raise ValueError(f"Variant position is outside reference sequence: {text!r}")
    if reference[position - 1] != wild_type:
        raise ValueError(
            f"Variant {token} has wild-type residue inconsistent with reference at position {position}"
        )
    if wild_type == mutant:
        raise ValueError(f"Synonymous variant is not a missense substitution: {text!r}")
    return wild_type, position, mutant


def _canonical_without_sequence(value: object) -> tuple[str, int, str] | None:
    if not isinstance(value, str):
        return None
    match = _VARIANT_RE.fullmatch(value.strip().upper())
    if match is None:
        return None
    wild_type, position_text, mutant = match.groups()
    if wild_type not in _AMINO_ACIDS or mutant not in _AMINO_ACIDS or wild_type == mutant:
        return None
    return wild_type, int(position_text), mutant


def _canonicalize(frame: pd.DataFrame, role: str) -> tuple[pd.DataFrame, int]:
    # Prefer an explicit full variant token.  Sequence-validated internal
    # frames also carry ``mutant`` as the one-letter replacement residue.
    column = "variant" if "variant" in frame.columns else "mutant" if "mutant" in frame.columns else None
    if column is None:
        raise ValueError(f"{role} assay is missing mutant/variant column")
    rows: list[dict[str, object]] = []
    for source_index, value in frame[column].items():
        try:
            key = _canonical_without_sequence(value)
        except (TypeError, ValueError):
            key = None
        if key is None:
            continue
        wt, position, mutant = key
        # Remove any pre-existing canonical fields before adding normalized
        # values; a raw ``mutant`` column often contains ``D2W`` while the
        # canonical ``mutant`` field must contain only ``W``.
        source_columns = [column for column in frame.columns if column not in {"wt", "position", "mutant"}]
        row = frame.loc[source_index, source_columns].to_dict()
        row.update({"wt": wt, "position": position, "mutant": mutant})
        rows.append(row)
    result = pd.DataFrame(rows)
    if result.empty:
        result = pd.DataFrame(columns=[*frame.columns, "wt", "position", "mutant"])
    # Callers may already have canonical columns (for example, a sequence-
    # validated assay table).  Keep the first occurrence rather than creating
    # duplicate labels that pandas cannot safely merge.
    result = result.loc[:, ~result.columns.duplicated()]
    if result.duplicated(["wt", "position", "mutant"]).any():
        duplicate_keys = result.loc[
            result.duplicated(["wt", "position", "mutant"], keep=False),
            ["wt", "position", "mutant"],
        ].drop_duplicates().to_dict("records")
        raise ValueError(f"{role} assay contains duplicate canonical variant keys: {duplicate_keys}")
    return result, len(result)


def pair_assays(abundance: pd.DataFrame, function: pd.DataFrame) -> pd.DataFrame:
    """Return an exact inner, one-to-one join of missense assay rows.

    Invalid/multi-mutant rows are excluded and recorded in ``DataFrame.attrs``.
    Duplicate canonical keys are fatal: averaging or a many-to-many merge would
    obscure whether a row is a replicate or an accidental duplicate.
    """
    if not isinstance(abundance, pd.DataFrame) or not isinstance(function, pd.DataFrame):
        raise TypeError("abundance and function must be pandas DataFrames")
    abundance_c, abundance_valid = _canonicalize(abundance, "abundance")
    function_c, function_valid = _canonicalize(function, "function")
    key = ["wt", "position", "mutant"]
    left = abundance_c.rename(columns={column: f"{column}_abundance" for column in abundance_c.columns if column not in key})
    right = function_c.rename(columns={column: f"{column}_function" for column in function_c.columns if column not in key})
    paired = left.merge(right, on=key, how="inner", validate="one_to_one", sort=True)
    abundance_keys = set(map(tuple, abundance_c[key].to_numpy().tolist()))
    function_keys = set(map(tuple, function_c[key].to_numpy().tolist()))
    paired.attrs["attrition"] = {
        "abundance_input": len(abundance),
        "function_input": len(function),
        "abundance_valid_single": abundance_valid,
        "function_valid_single": function_valid,
        "shared_pairs": len(paired),
        "abundance_unmatched": len(abundance_keys - function_keys),
        "function_unmatched": len(function_keys - abundance_keys),
    }
    return paired
