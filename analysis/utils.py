"""Shared helpers: loading universes, finding chain boundaries, and the
residue-index regions used to split the fusion construct into its chains
(Aga2-p, linker, RVG, nAChR).
"""
from __future__ import annotations

from pathlib import Path
from typing import Dict, List, Tuple, Union

import MDAnalysis as mda
import numpy as np
from MDAnalysis.core.groups import AtomGroup

PathLike = Union[str, Path]

THREE_TO_ONE: Dict[str, str] = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
    "HSD": "H", "HSE": "H", "HSP": "H",
}

# Residue-index boundaries for each construct, found by eyeballing sequences
# against the known design (Aga2-p -> linker -> RVG -> nAChR). The chains
# are numbered contiguously in these structures, so there's no resid gap to
# detect automatically - see detect_chain_boundaries below. None = to the end.
PHASE2_REGION_BOUNDS: Dict[str, Tuple[int, "int | None"]] = {
    "rvg": (0, 29),
    "nachr": (29, None),
}

PHASE3_REGION_BOUNDS: Dict[str, Tuple[int, "int | None"]] = {
    "aga2p": (0, 93),
    "linker": (93, 109),
    "rvg": (109, 139),
    "fusion": (0, 139),
    "nachr_a": (139, 461),
    "nachr_b": (506, 828),
    "nachr": (139, None),
}


def load_universe(gro_path: PathLike, xtc_path: PathLike) -> mda.Universe:
    """Load a structure + trajectory pair, with a clearer error if a path is wrong."""
    gro_path, xtc_path = Path(gro_path), Path(xtc_path)
    for p in (gro_path, xtc_path):
        if not p.exists():
            raise FileNotFoundError(f"Input file not found: {p}")
    try:
        return mda.Universe(str(gro_path), str(xtc_path))
    except Exception as exc:
        raise OSError(f"Failed to load universe from {gro_path}, {xtc_path}: {exc}") from exc


def detect_chain_boundaries(universe: mda.Universe) -> List[Tuple[int, int, int]]:
    """Find gaps in protein resid numbering - usually marks a chain break.

    Doesn't help here since our chains are numbered contiguously across the
    fusion construct; kept for structures where it would actually work.
    Returns (resid_before, resid_after, gap_size) tuples.
    """
    resids = universe.select_atoms("protein").residues.resids
    gaps = np.where(np.diff(resids) > 1)[0]
    return [
        (int(resids[g]), int(resids[g + 1]), int(resids[g + 1] - resids[g]))
        for g in gaps
    ]


def get_backbone_selection(atoms: AtomGroup) -> str:
    """Index-based MDAnalysis selection string for an atom group's backbone."""
    if len(atoms) == 0:
        raise ValueError("Cannot build a selection from an empty AtomGroup")
    return f"index {atoms.indices[0]}:{atoms.indices[-1]} and backbone"


def residue_sequence(residues) -> str:
    """One-letter sequence for a ResidueGroup; unknown resnames become '?'."""
    return "".join(THREE_TO_ONE.get(r.resname, "?") for r in residues)


def get_phase_regions(universe: mda.Universe, phase: int) -> Dict[str, AtomGroup]:
    """Split a universe's protein residues into named regions by index.

    phase=2 for RVG-only runs, phase=3 for the full fusion protein.
    """
    bounds = {2: PHASE2_REGION_BOUNDS, 3: PHASE3_REGION_BOUNDS}.get(phase)
    if bounds is None:
        raise ValueError(f"phase must be 2 or 3, got {phase!r}")

    residues = universe.select_atoms("protein").residues
    return {name: residues[start:stop].atoms for name, (start, stop) in bounds.items()}
