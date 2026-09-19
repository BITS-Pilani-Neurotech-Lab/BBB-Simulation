#!/usr/bin/env python
"""Place the receptor at a fixed distance from a ligand reference point and
combine both into a two-chain complex PDB.

The ligand (fusion protein) stays fixed; the receptor (nAChR) gets
translated along the ligand-to-binding-site vector, far enough apart to
avoid clashes but close enough for free MD to bring them together.
"""
from __future__ import annotations

import argparse
from pathlib import Path
from typing import Optional, Sequence

import numpy as np
from Bio.PDB import PDBIO, Chain, Model, PDBParser, Structure


def place_at_distance(
    ligand_path: Path,
    receptor_path: Path,
    output_path: Path,
    ligand_reference_point: np.ndarray,
    receptor_binding_center: np.ndarray,
    target_distance_angstrom: float = 30.0,
    binding_residues: Optional[Sequence[int]] = None,
) -> Path:
    for p in (ligand_path, receptor_path):
        if not p.exists():
            raise FileNotFoundError(f"Input file not found: {p}")

    parser = PDBParser(QUIET=True)
    ligand = parser.get_structure("ligand", str(ligand_path))
    receptor = parser.get_structure("receptor", str(receptor_path))

    direction = ligand_reference_point - receptor_binding_center
    current_distance = np.linalg.norm(direction)
    direction = direction / current_distance
    translation = direction * (current_distance - target_distance_angstrom)

    for atom in receptor.get_atoms():
        atom.coord = atom.coord + translation

    if binding_residues:
        # sanity check: distance from the given residues' CA to the reference point
        receptor_chain = list(receptor[0].get_chains())[0]
        ca_coords = [
            res["CA"].coord
            for res in receptor_chain.get_residues()
            if res.id[1] in binding_residues and "CA" in res
        ]
        if ca_coords:
            new_center = np.mean(ca_coords, axis=0)
            print(f"New distance: {np.linalg.norm(new_center - ligand_reference_point):.2f} A")

    combined = Structure.Structure("complex")
    model = Model.Model(0)
    chain_a, chain_b = Chain.Chain("A"), Chain.Chain("B")

    for res in list(ligand[0].get_chains())[0].get_residues():
        chain_a.add(res.copy())
    for res in list(receptor[0].get_chains())[0].get_residues():
        chain_b.add(res.copy())

    model.add(chain_a)
    model.add(chain_b)
    combined.add(model)

    output_path.parent.mkdir(parents=True, exist_ok=True)
    io = PDBIO()
    io.set_structure(combined)
    io.save(str(output_path))
    return output_path


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("ligand_pdb", type=Path)
    parser.add_argument("receptor_pdb", type=Path)
    parser.add_argument("output_pdb", type=Path)
    parser.add_argument("--ligand-point", type=float, nargs=3, required=True, metavar=("X", "Y", "Z"))
    parser.add_argument("--receptor-center", type=float, nargs=3, required=True, metavar=("X", "Y", "Z"))
    parser.add_argument("--distance", type=float, default=30.0)
    parser.add_argument("--binding-residues", type=int, nargs="*", default=None)
    args = parser.parse_args()

    out = place_at_distance(
        args.ligand_pdb,
        args.receptor_pdb,
        args.output_pdb,
        np.array(args.ligand_point),
        np.array(args.receptor_center),
        args.distance,
        args.binding_residues,
    )
    print(f"Saved {out}")


if __name__ == "__main__":
    main()
