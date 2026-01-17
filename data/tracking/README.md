# Ligand Tracking Pipeline

Sequential tracking system for GPR55 antagonist discovery.

## Pipeline Flow

```
01_unfiltered_library (35K+)
         ↓ Lipinski Filter
02_filtered_druglike (~20K)
         ↓ HTVS Docking
03_top10pct_per_target (P0, P3, Interface)
         ↓ Affinity Ranking
04_top250_per_target (250 × 3 targets)
         ↓ Selectivity + Drug Score
05_top25_per_target (25 × 3 targets)
         ↓ Expert Selection
06_top3_md_candidates (3 total)
         ↓ MD Validation
07_controls (AM251, ML193, etc.)
```

## Files

| Step | File | Description |
|------|------|-------------|
| 1 | `01_unfiltered_library.csv` | All 35K+ raw ligands from ChEMBL/PubChem |
| 2 | `02_filtered_druglike.csv` | Passed Lipinski Rule of 5 |
| 3 | `03_top10pct_per_target.csv` | Top 10% by affinity for each target |
| 4 | `04_top250_per_target.csv` | Top 250 per target (P0, P3, Interface) |
| 5 | `05_top25_per_target.csv` | Top 25 per target with selectivity scores |
| 6 | `06_top3_md_candidates.csv` | Final 3 for 100ns MD simulation |
| 7 | `07_controls.csv` | Reference antagonists (AM251, ML193) |

## Target Sites

| Site | Type | Description |
|------|------|-------------|
| P0 | Orthosteric | Primary binding pocket |
| P3 | Allosteric | PPI interface site |
| Interface | Hybrid | Membrane-receptor interface |

## Key Columns

- **affinity_kcal_mol**: AutoDock Vina docking score
- **selectivity_ratio**: P3/P0 ratio (>1 = P3 selective)
- **md_status**: pending / running / completed / failed
- **binding_free_energy**: MM-PBSA result from MD
