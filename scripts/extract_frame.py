#!/usr/bin/env python
"""Dump a single trajectory frame (default: last) to a standalone PDB."""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from analysis.utils import load_universe  # noqa: E402


def extract_frame(
    gro_path: Path,
    xtc_path: Path,
    output_path: Path,
    frame: int = -1,
    selection: str = "protein",
) -> Path:
    u = load_universe(gro_path, xtc_path)
    u.trajectory[frame]
    atoms = u.select_atoms(selection)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    atoms.write(str(output_path))
    return output_path


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("gro", type=Path)
    parser.add_argument("xtc", type=Path)
    parser.add_argument("output_pdb", type=Path)
    parser.add_argument("--frame", type=int, default=-1)
    parser.add_argument("--selection", default="protein")
    args = parser.parse_args()

    out = extract_frame(args.gro, args.xtc, args.output_pdb, args.frame, args.selection)
    print(f"Saved {out}")


if __name__ == "__main__":
    main()
