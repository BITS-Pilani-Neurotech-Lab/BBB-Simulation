#!/usr/bin/env python
"""Phase 2 pipeline: RVG-only MD analysis. Global + chain-resolved RMSD,
RVG-nAChR contact distances. Figures -> results/figures, arrays -> results/data.
"""
from __future__ import annotations

import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

import numpy as np  # noqa: E402

from analysis.contacts import compute_contact_distances, contact_fraction  # noqa: E402
from analysis.rmsd import compute_chain_resolved_rmsd, compute_global_rmsd, phase_decompose  # noqa: E402
from analysis.utils import get_phase_regions, load_universe  # noqa: E402
from analysis.visualization import (  # noqa: E402
    plot_chain_resolved_rmsd,
    plot_contact_distances,
    plot_global_rmsd,
)

GRO = PROJECT_ROOT / "data" / "trajectories" / "phase2" / "clean_stage.gro"
XTC = PROJECT_ROOT / "data" / "trajectories" / "phase2" / "final_smooth_movie.xtc"
FIGURES_DIR = PROJECT_ROOT / "results" / "figures"
DATA_DIR = PROJECT_ROOT / "results" / "data"


def main() -> None:
    u = load_universe(GRO, XTC)

    rmsd_data = compute_global_rmsd(u)
    np.save(DATA_DIR / "phase2_rmsd_global.npy", rmsd_data)
    time_ns = rmsd_data[:, 1] / 1000
    plot_global_rmsd(
        time_ns,
        rmsd_data[:, 2],
        FIGURES_DIR / "phase2_rmsd_trajectory.png",
        title="RVG-nAChR Backbone RMSD Over Simulation",
    )

    stats = phase_decompose(rmsd_data)
    for label, s in stats.items():
        print(f"{label}: {s['mean']:.3f} +/- {s['std']:.3f}")

    chain_rmsd = compute_chain_resolved_rmsd(u, phase=2)
    for name, arr in chain_rmsd.items():
        np.save(DATA_DIR / f"phase2_rmsd_{name}.npy", arr)
    plot_chain_resolved_rmsd(chain_rmsd, FIGURES_DIR / "phase2_rmsd_split.png")

    regions = get_phase_regions(u, phase=2)
    contact_data = compute_contact_distances(u, regions["rvg"], regions["nachr"])
    np.save(
        DATA_DIR / "phase2_contact_distances.npy",
        np.column_stack(
            [contact_data["time_ns"], contact_data["com_distance_nm"], contact_data["min_distance_nm"]]
        ),
    )
    plot_contact_distances(contact_data, FIGURES_DIR / "phase2_contact_analysis.png")

    frac = contact_fraction(contact_data["min_distance_nm"])
    print(
        f"CoM distance: start={contact_data['com_distance_nm'][0]:.3f} "
        f"end={contact_data['com_distance_nm'][-1]:.3f} nm"
    )
    print(
        f"Min distance: start={contact_data['min_distance_nm'][0]:.3f} "
        f"end={contact_data['min_distance_nm'][-1]:.3f} nm"
    )
    print(f"Fraction of frames in contact (<0.35 nm): {frac:.1%}")


if __name__ == "__main__":
    main()
