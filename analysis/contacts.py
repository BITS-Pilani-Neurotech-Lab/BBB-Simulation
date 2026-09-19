"""Contact-distance analysis: CoM/min-distance over a trajectory, plus
per-residue contacts on a single frame.
"""
from __future__ import annotations

from typing import Dict

import MDAnalysis as mda
import numpy as np
import pandas as pd
from MDAnalysis.analysis import distances
from MDAnalysis.core.groups import AtomGroup, ResidueGroup
from MDAnalysis.lib.distances import calc_bonds

ANGSTROM_TO_NM = 0.1


def compute_contact_distances(
    universe: mda.Universe, region_a: AtomGroup, region_b: AtomGroup
) -> Dict[str, np.ndarray]:
    """Per-frame CoM distance and minimum heavy-atom distance between two regions.

    Walks the whole trajectory, so position/align `universe` beforehand if needed.
    Distances use the minimum-image convention (box-aware) so a molecule
    getting PBC-wrapped to a different periodic image on some frame - which
    `md_centered.xtc` does on at least one frame - doesn't show up as a
    spurious multi-nm spike. Returns time_ns / com_distance_nm / min_distance_nm.
    """
    times, com_distances, min_distances = [], [], []
    for ts in universe.trajectory:
        times.append(ts.time / 1000)
        box = ts.dimensions
        com_a = region_a.center_of_mass().reshape(1, 3)
        com_b = region_b.center_of_mass().reshape(1, 3)
        com_distances.append(calc_bonds(com_a, com_b, box=box)[0] * ANGSTROM_TO_NM)
        dist_matrix = distances.distance_array(region_a.positions, region_b.positions, box=box)
        min_distances.append(dist_matrix.min() * ANGSTROM_TO_NM)

    return {
        "time_ns": np.array(times),
        "com_distance_nm": np.array(com_distances),
        "min_distance_nm": np.array(min_distances),
    }


def contact_fraction(min_distances: np.ndarray, threshold: float = 0.35) -> float:
    """Fraction of frames with min distance below `threshold` (nm)."""
    return float((np.asarray(min_distances) < threshold).mean())


def compute_per_residue_contacts(
    region: AtomGroup,
    region_residues: ResidueGroup,
    reference: AtomGroup,
    threshold: float = 0.35,
) -> pd.DataFrame:
    """Which residues of `region` have an atom within `threshold` nm of `reference`.

    Uses current atom positions, so this is normally run on one frame (e.g.
    after seeking to the last frame). Box-aware (minimum image), same reason
    as compute_contact_distances. Returns resname/resid/min_distance_nm/in_contact.
    """
    box = region.universe.dimensions
    dist_matrix = distances.distance_array(region.positions, reference.positions, box=box) * ANGSTROM_TO_NM
    base_index = region.indices[0]

    rows = []
    for res in region_residues:
        res_rows = [a.index - base_index for a in res.atoms if a.index in region.indices]
        if not res_rows:
            continue
        min_dist = float(dist_matrix[res_rows, :].min())
        rows.append(
            {
                "resname": res.resname,
                "resid": int(res.resid),
                "min_distance_nm": min_dist,
                "in_contact": min_dist < threshold,
            }
        )
    return pd.DataFrame(rows)
