# Molecular Dynamics Simulations

This folder contains MD simulation inputs and analysis for top HTVS hits.

## Structure
```
simulations/
├── inputs/           # Topology, coordinates, .mdp parameter files
├── analysis/         # RMSD, RMSF, MMPBSA results (.xvg, .csv)
└── README.md
```

## Files Tracked by Git
- Input files: `*.gro`, `*.top`, `*.mdp`, `*.pdb`
- Analysis outputs: `*_analysis.xvg`, `*.csv`

## Files NOT Tracked (too large)
Trajectory and checkpoint files are excluded via `.gitignore`:
- `*.xtc`, `*.trr`, `*.dcd` (trajectory)
- `*.cpt`, `*.rst7` (checkpoints)
- `*.edr` (energy)

> **Note**: Store trajectories in Google Drive, GCP Storage, or Zenodo for archival.

## Simulation Protocol
- Software: GROMACS 2023+
- Force Field: CHARMM36m
- Duration: 100 ns per system
- Ensemble: NPT (300K, 1 bar)
