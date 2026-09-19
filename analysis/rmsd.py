"""RMSD analysis - global and per-chain backbone RMSD, alignment, phase stats."""
from __future__ import annotations

from typing import Dict

import MDAnalysis as mda
import numpy as np
from MDAnalysis.analysis import align, rms

from .utils import get_backbone_selection, get_phase_regions


def compute_global_rmsd(
    universe: mda.Universe, select: str = "backbone", ref_frame: int = 0
) -> np.ndarray:
    """RMSD of `select` against `ref_frame`. Returns [frame, time_ps, rmsd_angstrom] columns."""
    R = rms.RMSD(universe, select=select, ref_frame=ref_frame)
    R.run()
    return R.results.rmsd


def compute_chain_resolved_rmsd(
    universe: mda.Universe, phase: int, ref_frame: int = 0
) -> Dict[str, np.ndarray]:
    """Backbone RMSD per region (see utils.get_phase_regions) instead of the whole thing."""
    regions = get_phase_regions(universe, phase)
    results: Dict[str, np.ndarray] = {}
    for name, atoms in regions.items():
        R = rms.RMSD(universe, select=get_backbone_selection(atoms), ref_frame=ref_frame)
        R.run()
        results[name] = R.results.rmsd
    return results


def align_trajectory(
    universe: mda.Universe, reference_selection: str, in_memory: bool = True
) -> mda.Universe:
    """Align the trajectory to `reference_selection`'s first frame, in place.

    in_memory=True is needed for the alignment to stick around across later
    analysis passes over the same universe.
    """
    aligner = align.AlignTraj(universe, universe, select=reference_selection, in_memory=in_memory)
    aligner.run()
    return universe


def phase_decompose(rmsd_array: np.ndarray, n_quarters: int = 4) -> Dict[str, Dict[str, float]]:
    """Mean +/- std over the first/second half and n_quarters equal chunks of an RMSD trace."""
    vals = rmsd_array[:, 2] if rmsd_array.ndim == 2 else np.asarray(rmsd_array)
    n = len(vals)
    if n == 0:
        raise ValueError("rmsd_array is empty")

    midpoint = n // 2
    stats: Dict[str, Dict[str, float]] = {
        "first_half": {"mean": float(vals[:midpoint].mean()), "std": float(vals[:midpoint].std())},
        "second_half": {"mean": float(vals[midpoint:].mean()), "std": float(vals[midpoint:].std())},
    }

    q = n // n_quarters
    for i in range(n_quarters):
        start = i * q
        stop = (i + 1) * q if i < n_quarters - 1 else n
        seg = vals[start:stop]
        stats[f"q{i + 1}"] = {"mean": float(seg.mean()), "std": float(seg.std())}
    return stats
