#!/usr/bin/env python
"""Phase 3 pipeline: full fusion-protein MD analysis. Computes RVG-nAChR
contact distances and per-residue contacts on the raw trajectory, then
aligns to the more stable nAChR chain for chain-resolved RMSD. Figures ->
results/figures, arrays/CSVs -> results/data.

Contact/distance work has to happen before alignment: AlignTraj rotates
atom coordinates into a best-fit frame but leaves the box vectors as-is, so
periodic (minimum-image) distance calculations done afterward are wrong.
"""
from __future__ import annotations

import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

import numpy as np  # noqa: E402

from analysis.contacts import compute_contact_distances, compute_per_residue_contacts, contact_fraction  # noqa: E402
from analysis.rmsd import align_trajectory, compute_chain_resolved_rmsd  # noqa: E402
from analysis.utils import PHASE3_REGION_BOUNDS, get_backbone_selection, get_phase_regions, load_universe  # noqa: E402
from analysis.visualization import plot_chain_resolved_rmsd, plot_contact_distances  # noqa: E402

GRO = PROJECT_ROOT / "data" / "trajectories" / "phase3" / "md_fusion.gro"
XTC = PROJECT_ROOT / "data" / "trajectories" / "phase3" / "md_centered.xtc"
FIGURES_DIR = PROJECT_ROOT / "results" / "figures"
DATA_DIR = PROJECT_ROOT / "results" / "data"
CONTACT_THRESHOLD_NM = 0.35


def main() -> None:
    u = load_universe(GRO, XTC)
    regions = get_phase_regions(u, phase=3)
    nachr_atoms = regions["nachr"]

    contact_data = compute_contact_distances(u, regions["rvg"], nachr_atoms)
    np.save(
        DATA_DIR / "phase3_contact_distances.npy",
        np.column_stack(
            [contact_data["time_ns"], contact_data["com_distance_nm"], contact_data["min_distance_nm"]]
        ),
    )
    plot_contact_distances(
        contact_data, FIGURES_DIR / "phase3_contact_analysis.png", title_suffix=" (Full Fusion Simulation)"
    )
    frac = contact_fraction(contact_data["min_distance_nm"], CONTACT_THRESHOLD_NM)
    print(f"Fraction of frames in contact (<{CONTACT_THRESHOLD_NM} nm): {frac:.1%}")

    # per-residue contacts at the final frame only
    u.trajectory[-1]
    all_residues = u.select_atoms("protein").residues
    for label in ("aga2p", "rvg"):
        start, stop = PHASE3_REGION_BOUNDS[label]
        region_atoms = regions[label]
        region_residues = all_residues[start:stop]
        df = compute_per_residue_contacts(region_atoms, region_residues, nachr_atoms, CONTACT_THRESHOLD_NM)
        df.to_csv(DATA_DIR / f"phase3_per_residue_contacts_{label}.csv", index=False)
        n_contact = int(df["in_contact"].sum())
        print(f"{label}: {n_contact}/{len(df)} residues in contact with nAChR")

    # RMSD needs a common reference frame, so it's the one thing that runs post-alignment
    align_trajectory(u, get_backbone_selection(regions["nachr_b"]))

    chain_rmsd = compute_chain_resolved_rmsd(u, phase=3)
    plotted_regions = ("fusion", "rvg", "nachr_a", "nachr_b")
    for name in plotted_regions:
        np.save(DATA_DIR / f"phase3_rmsd_{name}.npy", chain_rmsd[name])
    plot_chain_resolved_rmsd(
        {name: chain_rmsd[name] for name in plotted_regions},
        FIGURES_DIR / "phase3_rmsd_aligned.png",
    )


if __name__ == "__main__":
    main()
