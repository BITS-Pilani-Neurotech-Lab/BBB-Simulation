#!/usr/bin/env python
"""Strip a PDB down to the record types HADDOCK accepts."""
from __future__ import annotations

import argparse
from pathlib import Path
from typing import Tuple

VALID_RECORDS = ("ATOM", "HETATM", "TER", "END", "REMARK", "HEADER", "SEQRES", "SSBOND", "CRYST1")


def clean_pdb_for_haddock(input_path: Path, output_path: Path) -> Tuple[int, int]:
    if not input_path.exists():
        raise FileNotFoundError(f"Input file not found: {input_path}")

    with open(input_path) as f:
        lines = f.readlines()

    cleaned = [line for line in lines if line.strip() and line[:6].strip() in VALID_RECORDS]
    if not cleaned:
        raise ValueError(f"No valid PDB records found in {input_path}")
    if not cleaned[-1].startswith("END"):
        cleaned.append("END\n")

    output_path.parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, "w") as f:
        f.writelines(cleaned)

    return len(lines), len(cleaned)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input_pdb", type=Path)
    parser.add_argument("output_pdb", type=Path)
    args = parser.parse_args()

    n_in, n_out = clean_pdb_for_haddock(args.input_pdb, args.output_pdb)
    print(f"{args.input_pdb}: {n_in} lines -> {args.output_pdb}: {n_out} lines")


if __name__ == "__main__":
    main()
