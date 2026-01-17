# Ligand Tracking System

Comprehensive tracking sheets for the GPR55 antagonist discovery campaign.

## Files

| Sheet | Description | When to Update |
|-------|-------------|----------------|
| `01_all_ligands.csv` | All 35K+ ligands with IDs, SMILES, sources | After library preparation |
| `02_filtered_ligands.csv` | Lipinski/drug-like filtered compounds | After filtering step |
| `03_htvs_results.csv` | Docking results (affinity, rank, hit category) | After each HTVS batch |
| `04_top25_hits.csv` | Top 25 ranked compounds across all targets | After HTVS completion |
| `05_md_candidates.csv` | Top 3 for MD simulation with tracking | During/after MD runs |

## Column Definitions

### Common Fields
| Column | Description |
|--------|-------------|
| `ligand_id` | Unique identifier (CHEMBL ID or PubChem CID) |
| `source` | Database: ChEMBL, PubChem, ZINC, etc. |
| `source_id` | ID in original database |
| `smiles` | Canonical SMILES string |
| `batch` | Batch number (01-09) |

### Molecular Properties
| Column | Description |
|--------|-------------|
| `molecular_weight` | MW in Da |
| `hbd` | H-bond donors |
| `hba` | H-bond acceptors |
| `logp` | Partition coefficient |
| `tpsa` | Topological polar surface area |
| `rotatable_bonds` | Number of rotatable bonds |
| `lipinski_pass` | True/False for Lipinski Rule of 5 |

### Docking Results
| Column | Description |
|--------|-------------|
| `target_site` | P0 (Orthosteric), P3 (Allosteric), Interface |
| `affinity_kcal_mol` | Best docking score |
| `percentile` | top_1%, top_5%, top_10%, etc. |
| `hit_category` | super_hit (<-10), hit (<-8.5), moderate (<-7) |
| `selectivity_p3_p0` | Ratio for allosteric selectivity |

### MD Simulation
| Column | Description |
|--------|-------------|
| `md_status` | pending, running, completed, failed |
| `md_duration_ns` | Simulation length in nanoseconds |
| `rmsd_avg` | Average RMSD (Å) |
| `rmsf_binding_site` | RMSF at binding site (Å) |
| `binding_free_energy_kcal_mol` | MM-PBSA result |

## Workflow

```
01_all_ligands → 02_filtered_ligands → 03_htvs_results
                                              ↓
                                       04_top25_hits
                                              ↓
                                       05_md_candidates
```
