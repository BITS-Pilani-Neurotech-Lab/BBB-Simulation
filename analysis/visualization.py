"""Plotting for RMSD and contact-distance results."""
from __future__ import annotations

from pathlib import Path
from typing import Dict, Optional, Union

import matplotlib.pyplot as plt
import numpy as np

PathLike = Union[str, Path]


def _save(fig: plt.Figure, output_path: PathLike) -> None:
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=300)


def plot_global_rmsd(
    time_ns: np.ndarray,
    rmsd_vals: np.ndarray,
    output_path: PathLike,
    title: str = "Backbone RMSD Over Simulation",
) -> plt.Figure:
    """Single RMSD trace with equilibration/production halves shaded."""
    midpoint = len(time_ns) // 2
    production_mean = rmsd_vals[midpoint:].mean()

    fig, ax = plt.subplots(figsize=(10, 4))
    ax.plot(time_ns, rmsd_vals, color="steelblue", linewidth=0.8, alpha=0.9)
    ax.axvspan(time_ns[0], time_ns[midpoint], alpha=0.1, color="red", label="Equilibration phase")
    ax.axvspan(time_ns[midpoint], time_ns[-1], alpha=0.1, color="green", label="Production phase")
    ax.axhline(
        production_mean,
        color="green",
        linestyle="--",
        linewidth=1.2,
        label=f"Production mean: {production_mean:.2f}",
    )
    ax.set_xlabel("Time (ns)")
    ax.set_ylabel("Backbone RMSD")
    ax.set_title(title)
    ax.legend()
    fig.tight_layout()
    _save(fig, output_path)
    return fig


def plot_chain_resolved_rmsd(
    rmsd_by_region: Dict[str, np.ndarray],
    output_path: PathLike,
    colors: Optional[Dict[str, str]] = None,
) -> plt.Figure:
    """One RMSD panel per region, stacked vertically. `rmsd_by_region` from rmsd.compute_chain_resolved_rmsd."""
    default_colors = ["steelblue", "darkorange", "green", "purple", "brown", "teal"]
    colors = colors or {}

    names = list(rmsd_by_region.keys())
    fig, axes = plt.subplots(len(names), 1, figsize=(10, 3 * len(names)), sharex=True)
    axes = np.atleast_1d(axes)

    time_ns = next(iter(rmsd_by_region.values()))[:, 1] / 1000
    midpoint = len(time_ns) // 2

    for i, (ax, name) in enumerate(zip(axes, names)):
        color = colors.get(name, default_colors[i % len(default_colors)])
        vals = rmsd_by_region[name][:, 2]
        mean = vals[midpoint:].mean()
        ax.plot(time_ns, vals, color=color, linewidth=0.8)
        ax.axhline(
            mean, color="red", linestyle="--", linewidth=1.2, label=f"Production mean: {mean:.3f} nm"
        )
        ax.set_ylabel("RMSD (nm)")
        ax.set_title(name)
        ax.legend(fontsize=9)

    axes[-1].set_xlabel("Time (ns)")
    fig.tight_layout()
    _save(fig, output_path)
    return fig


def plot_contact_distances(
    contact_data: Dict[str, np.ndarray],
    output_path: PathLike,
    contact_threshold: float = 0.35,
    title_suffix: str = "",
) -> plt.Figure:
    """CoM distance + min heavy-atom distance panels. `contact_data` from contacts.compute_contact_distances."""
    times = contact_data["time_ns"]
    com_distances = contact_data["com_distance_nm"]
    min_distances = contact_data["min_distance_nm"]

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 7), sharex=True)

    ax1.plot(times, com_distances, color="purple", linewidth=0.8)
    ax1.axhline(2.0, color="red", linestyle="--", linewidth=1, label="2 nm reference")
    ax1.set_ylabel("CoM Distance (nm)")
    ax1.set_title(f"RVG-nAChR Center of Mass Distance{title_suffix}")
    ax1.legend()

    ax2.plot(times, min_distances, color="green", linewidth=0.8)
    ax2.axhline(
        contact_threshold,
        color="red",
        linestyle="--",
        linewidth=1,
        label=f"Contact threshold ({contact_threshold} nm)",
    )
    ax2.set_ylabel("Min Distance (nm)")
    ax2.set_xlabel("Time (ns)")
    ax2.set_title(f"RVG-nAChR Minimum Heavy Atom Distance{title_suffix}")
    ax2.legend()

    fig.tight_layout()
    _save(fig, output_path)
    return fig
