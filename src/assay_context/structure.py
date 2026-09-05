"""Reproducible residue-level structural mapping and feature extraction.

The parser intentionally supports the small, text PDB format without making
Biopython a runtime requirement.  Experimental confidence is retained as the
mean atom B-factor; it is never renamed to pLDDT.  Missing residues remain
missing throughout the feature table.
"""

from __future__ import annotations

import math
import re
from collections.abc import Iterable
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd

AA3 = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
}
VDW = {"H": 1.20, "C": 1.70, "N": 1.55, "O": 1.52, "S": 1.80, "P": 1.80}
# Tien et al. 2013 residue maximum ASA values (square Angstroms).
MAX_ASA = {
    "A": 129, "R": 274, "N": 195, "D": 193, "C": 167, "Q": 225, "E": 223,
    "G": 104, "H": 224, "I": 197, "L": 201, "K": 236, "M": 224, "F": 240,
    "P": 159, "S": 155, "T": 172, "W": 285, "Y": 263, "V": 174,
}


@dataclass
class _Atom:
    name: str
    element: str
    xyz: np.ndarray
    bfactor: float


@dataclass
class _Residue:
    key: tuple[int, str]
    name: str
    aa: str
    atoms: list[_Atom]


def _parse_pdb(path: str, chain: str) -> tuple[list[_Residue], dict[tuple[int, str], str], str | None]:
    residues: dict[tuple[int, str], _Residue] = {}
    helix: dict[tuple[int, str], str] = {}
    found_chain = False
    active_model = 1
    model = 1
    for raw in Path(path).read_text(encoding="utf-8", errors="replace").splitlines():
        record = raw[:6].strip().upper()
        if record == "MODEL":
            try:
                model = int(raw[10:14].strip() or 1)
            except ValueError:
                model = 1
            continue
        if record == "ENDMDL":
            if model == active_model:
                break
            continue
        if model != active_model or record not in {"ATOM", "HELIX", "SHEET"}:
            continue
        if record in {"HELIX", "SHEET"}:
            # HELIX/SHEET records use fixed-width chain and residue fields.
            if record == "HELIX":
                c1, n1, c2, n2 = raw[19].strip(), raw[21:25].strip(), raw[31].strip(), raw[33:37].strip()
                code = "H"
            else:
                c1, n1, c2, n2 = raw[21].strip(), raw[22:26].strip(), raw[32].strip(), raw[33:37].strip()
                code = "E"
            if c1 == chain and c2 == chain:
                try:
                    for number in range(int(n1), int(n2) + 1):
                        helix[(number, "")] = code
                except ValueError:
                    pass
            continue
        c = raw[21].strip() or " "
        if c != chain:
            continue
        found_chain = True
        altloc = raw[16].strip()
        if altloc not in {"", "A", "1"}:
            continue
        name = raw[17:20].strip().upper()
        aa = AA3.get(name)
        if aa is None:
            continue
        try:
            number = int(raw[22:26].strip())
            xyz = np.array([float(raw[30:38]), float(raw[38:46]), float(raw[46:54])], dtype=float)
            bfactor = float(raw[60:66].strip() or "nan")
        except (ValueError, IndexError):
            continue
        insertion = raw[26].strip()
        key = (number, insertion)
        residue = residues.setdefault(key, _Residue(key, name, aa, []))
        element = raw[76:78].strip().upper() or re.sub(r"[^A-Z]", "", raw[12:16]).upper()[:1]
        if element == "H":
            continue
        residue.atoms.append(_Atom(raw[12:16].strip(), element, xyz, bfactor))
    if not found_chain:
        raise ValueError(f"requested chain {chain!r} was not found in structure")
    ordered = [residue for residue in sorted(residues.values(), key=lambda r: r.key) if residue.atoms]
    if not ordered:
        raise ValueError(f"requested chain {chain!r} contains no standard observed residues")
    return ordered, helix, chain


def _align(reference: str, observed: str) -> tuple[list[int | None], float, float]:
    """Global Needleman-Wunsch alignment, returning observed index per ref pos."""
    n, m = len(reference), len(observed)
    gap, match, mismatch = -1.0, 2.0, -1.0
    scores = np.full((n + 1, m + 1), -np.inf)
    scores[0, :] = np.arange(m + 1) * gap
    scores[:, 0] = np.arange(n + 1) * gap
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            scores[i, j] = max(
                scores[i - 1, j - 1] + (match if reference[i - 1] == observed[j - 1] else mismatch),
                scores[i - 1, j] + gap,
                scores[i, j - 1] + gap,
            )
    mapping: list[int | None] = [None] * n
    i, j = n, m
    while i or j:
        diagonal = i and j and np.isclose(scores[i, j], scores[i - 1, j - 1] + (match if reference[i - 1] == observed[j - 1] else mismatch))
        # Prefer a gap in the observed sequence on ties.  This deterministically
        # maps the final matching residue in A-G-G/A-G to reference position 3.
        up = i and np.isclose(scores[i, j], scores[i - 1, j] + gap)
        if diagonal:
            mapping[i - 1] = j - 1
            i -= 1
            j -= 1
        elif up:
            i -= 1
        else:
            j -= 1
    aligned = [idx for idx in mapping if idx is not None]
    identity = sum(reference[pos] == observed[idx] for pos, idx in enumerate(mapping) if idx is not None) / max(len(aligned), 1)
    coverage = len(aligned) / max(n, 1)
    return mapping, float(identity), float(coverage)


def _sphere_points(n: int = 96) -> np.ndarray:
    i = np.arange(n, dtype=float)
    z = 1.0 - 2.0 * (i + 0.5) / n
    radius = np.sqrt(np.maximum(0.0, 1.0 - z * z))
    theta = math.pi * (3.0 - math.sqrt(5.0)) * i
    return np.column_stack((radius * np.cos(theta), radius * np.sin(theta), z))


def _relative_sasa(residue: _Residue, all_atoms: list[_Atom], points: np.ndarray) -> float:
    accessible_counts: list[int] = []
    atoms = residue.atoms
    all_coords = np.array([atom.xyz for atom in all_atoms])
    for atom in atoms:
        radius = VDW.get(atom.element, 1.7) + 1.4
        surface = atom.xyz + points * radius
        blocked = np.zeros(len(surface), dtype=bool)
        # First select atoms that could occlude a point.  This keeps the
        # deterministic Shrake-Rupley calculation tractable for a 3--4k atom
        # chain while retaining the exact distance rule.
        distances = np.linalg.norm(all_coords - atom.xyz, axis=1)
        candidates = np.flatnonzero(distances < radius + 3.5)
        for index in candidates:
            other = all_atoms[index]
            if other is atom:
                continue
            other_radius = VDW.get(other.element, 1.7) + 1.4
            blocked |= np.einsum("ij,ij->i", surface - other.xyz, surface - other.xyz) < other_radius**2
            if blocked.all():
                break
        accessible_counts.append(int((~blocked).sum()))
    area = sum(
        4.0 * math.pi * (VDW.get(atom.element, 1.7) + 1.4) ** 2 * count / max(len(points), 1)
        for atom, count in zip(atoms, accessible_counts, strict=True)
    )
    return float(area / MAX_ASA[residue.aa])


def map_structure(
    sequence: str,
    structure_path: str,
    chain: str,
    functional_positions: Iterable[int] | None = None,
    functional_site_source: str | None = None,
    structure_kind: str | None = None,
) -> pd.DataFrame:
    """Map every reference position to a PDB chain and calculate features.

    Mapping accepts only chains with at least 95% aligned identity and 65%
    reference coverage.  Functional positions are supplied by an independent,
    predeclared annotation; proximity is descriptive and does not imply
    allostery.  Experimental structures use B-factor confidence and retain a
    NaN ``plddt`` column by design.
    """
    if not isinstance(sequence, str) or not sequence.strip():
        raise ValueError("reference sequence must be non-empty")
    residues, secondary, _ = _parse_pdb(structure_path, chain)
    reference = sequence.strip().upper()
    observed = "".join(residue.aa for residue in residues)
    mapping, identity, coverage = _align(reference, observed)
    if identity < 0.95:
        raise ValueError(f"structure chain identity {identity:.3f} is below 0.95")
    if coverage < 0.65:
        raise ValueError(f"structure chain coverage {coverage:.3f} is below 0.65")
    stem = Path(structure_path).stem
    is_alphafold = structure_kind == "alphafold" or stem.upper().startswith("AF-")
    source_match = re.search(r"(?:^|[_-])(\d[0-9A-Za-z]{3})(?:$|[_-])", stem, re.IGNORECASE)
    source = (f"AlphaFold DB model {stem}" if is_alphafold else f"experimental PDB {source_match.group(1).upper()}" if source_match else "experimental PDB")
    points = _sphere_points()
    all_atoms = [atom for residue in residues for atom in residue.atoms]
    functional = {int(p) for p in (functional_positions or ())}
    mapped_site_atoms = [residues[idx].atoms for pos, idx in enumerate(mapping, start=1) if idx is not None and pos in functional]
    flat_site_atoms = [atom for atoms in mapped_site_atoms for atom in atoms]
    observed_to_reference = {observed_index: position for position, observed_index in enumerate(mapping, start=1) if observed_index is not None}
    rows: list[dict[str, object]] = []
    for position, observed_index in enumerate(mapping, start=1):
        base: dict[str, object] = {
            "position": position, "reference_residue": reference[position - 1], "structure_source": source,
            "mapping_identity": identity, "mapping_coverage": coverage, "chain": chain,
            "functional_site_source": functional_site_source,
            "plddt": np.nan, "confidence_type": "pLDDT" if is_alphafold else "B-factor", "confidence_value": np.nan,
            "accessibility": np.nan, "relative_sasa": np.nan, "secondary_structure": pd.NA,
            "contacts": np.nan, "site_distance": np.nan, "mapping_status": "missing",
        }
        if observed_index is not None:
            residue = residues[observed_index]
            base.update({
                "mapping_status": "mapped", "observed_residue": residue.aa,
                "observed_resseq": residue.key[0], "insertion_code": residue.key[1],
                "accessibility": _relative_sasa(residue, all_atoms, points),
                "secondary_structure": next((code for key, code in secondary.items() if key[0] == residue.key[0]), "C"),
                "confidence_value": float(np.nanmean([atom.bfactor for atom in residue.atoms])),
            })
            if is_alphafold:
                base["plddt"] = base["confidence_value"]
            neighbors = []
            for other_index, other in enumerate(residues):
                other_position = observed_to_reference.get(other_index)
                if other_index == observed_index or (other_position is not None and abs(other_position - position) <= 1):
                    continue
                distances = np.linalg.norm(np.array([a.xyz for a in residue.atoms])[:, None] - np.array([a.xyz for a in other.atoms])[None, :], axis=2)
                if np.nanmin(distances) <= 8.0:
                    neighbors.append(other_index)
            base["contacts"] = float(len(neighbors))
            if flat_site_atoms and position not in functional:
                distances = np.linalg.norm(np.array([a.xyz for a in residue.atoms])[:, None] - np.array([a.xyz for a in flat_site_atoms])[None, :], axis=2)
                base["site_distance"] = float(np.nanmin(distances))
            elif position in functional:
                base["site_distance"] = 0.0
        rows.append(base)
    result = pd.DataFrame(rows)
    result.attrs.update({"mapping_identity": identity, "mapping_coverage": coverage, "structure_source": source, "functional_site_source": functional_site_source, "confidence_definition": "AlphaFold pLDDT from atom B-factor" if is_alphafold else "experimental atom B-factor mean; pLDDT is NaN"})
    return result
