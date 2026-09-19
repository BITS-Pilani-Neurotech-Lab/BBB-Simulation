# BBB Simulation: RVG-nAChR Binding at the Blood-Brain Barrier

MD analysis of a Rabies Virus Glycoprotein (RVG) peptide, and later a full
Aga2-p/linker/RVG fusion protein, binding to the nicotinic acetylcholine
receptor (nAChR) as a candidate route across the blood-brain barrier.

## Phases

1. **Phase 1 — fusion protein design.** AlphaFold/ColabFold model, cleaned
   up in PyMOL, validated with MolProbity. No simulation code here beyond
   `scripts/clean_pdb.py` (HADDOCK record cleanup) and
   `scripts/prepare_complex.py` (initial docking placement).
2. **Phase 2 — RVG-only MD.** `data/trajectories/phase2/`, run via
   `scripts/run_phase2.py`.
3. **Phase 3 — full fusion protein MD.** `data/trajectories/phase3/`, run
   via `scripts/run_phase3.py`.

## Layout

```
BBB Sim/
├── requirements.txt
├── data/
│   ├── structures/     fusion.pdb, fusion_haddock.pdb, complex, receptor, last frame
│   ├── trajectories/   phase2/, phase3/ .gro + .xtc
│   └── molprobity/     MolProbity validation reports (PDFs)
├── analysis/
│   ├── utils.py          universe loading, chain regions, sequence helpers
│   ├── rmsd.py            global + chain-resolved RMSD, alignment, phase stats
│   ├── contacts.py       CoM/min-distance + per-residue contact analysis
│   └── visualization.py   plotting
├── scripts/
│   ├── clean_pdb.py
│   ├── prepare_complex.py
│   ├── extract_frame.py
│   ├── run_phase2.py
│   └── run_phase3.py
├── results/
│   ├── figures/     RMSD/contact plots
│   └── data/         .npy / .csv outputs
└── gromacs/
    ├── mdp/           ions/minim/nvt/npt/md
    ├── topol.top, index.ndx
    └── README.md      gmx command reference
```

## Structures

`fusion.pdb` is the AlphaFold model as cleaned in PyMOL (chain A assigned,
protonated with `reduce`) — the structure everything else is built from.
`fusion_haddock.pdb` is that same structure stripped to HADDOCK-valid
records; it's what `prepare_complex.py` loads as the ligand when building
`complex_30A.pdb`. `protein2.pdb` is the nAChR receptor. `last_clean.pdb` /
`last_frame_protein.pdb` are the final MD frame.

`data/molprobity/` holds MolProbity validation PDFs for the structures
above (currently just the fusion model's Ramachandran report) — drop more
in there as you generate them.

## Setup

```bash
python -m venv venv
venv\Scripts\activate
pip install -r requirements.txt
```

## Running

```bash
python scripts/run_phase2.py
python scripts/run_phase3.py
```

Figures land in `results/figures/`, arrays/CSVs in `results/data/`; RMSD
phase means, contact fractions, and per-residue contact counts print to
stdout.

## Construct regions

RVG-only and fusion constructs are split into named regions by **residue
index**, not a detected chain gap — the chains are numbered contiguously,
so `detect_chain_boundaries` won't find anything. Boundaries were worked
out by eyeballing sequences against the known construct design
(`PHASE2_REGION_BOUNDS` / `PHASE3_REGION_BOUNDS` in `analysis/utils.py`):

- **Phase 2**: `rvg` = 0–28, `nachr` = 29–end.
- **Phase 3**: `aga2p` = 0–92, `linker` = 93–108, `rvg` = 109–138,
  `nachr_a` = 139–460, `nachr_b` = 506–827.

A structure with different residue ordering needs these re-derived.

## PBC gotcha

Contact-distance calculations (`analysis/contacts.py`) use minimum-image
distances (box-aware), since `md_centered.xtc` has at least one frame where
a region gets wrapped to a different periodic image. That only works if the
box is still in its original orientation — `align_trajectory`'s rigid-body
fit rotates atom coordinates but doesn't rotate the box vectors, so running
contacts on an already-aligned universe silently reintroduces the same
jump. `run_phase3.py` computes contacts and per-residue contacts *before*
aligning for RMSD, on purpose - keep that order if you touch it.

## Other notes

- No HADDOCK `.tbl`/`.cns` restraint files exist for this project.
- `last_frame_clean.pdb` and `last_protein.pdb` (dropped during cleanup)
  were byte-identical to `last_clean.pdb` and `last_frame_protein.pdb`.
  `complex_manual.pdb` and `complex_30A.pdb` are *not* identical despite
  matching sizes — both are kept.
- `gromacs/mdp/`, `topol.top`, and `index.ndx` came from the WSL GROMACS
  run directory.
