# GPR55 Allosteric Antagonist Discovery

> **Computational Drug Discovery Pipeline for Novel GPR55 Antagonists**

[![Phase](https://img.shields.io/badge/Phase-2%20HTVS-blue)](docs/PHASE_2_PREPARATION.md)
[![Stage](https://img.shields.io/badge/Stage-2%20Running-orange)]()
[![Target](https://img.shields.io/badge/Target-GPR55-green)](docs/GPR55_INTRODUCTION.md)
[![License](https://img.shields.io/badge/License-MIT-yellow)](LICENSE)

---

## 🎯 Project Overview

This project aims to discover novel allosteric antagonists for **GPR55** (G protein-coupled receptor 55), an atypical cannabinoid receptor implicated in cancer progression, inflammation, and metabolic disorders.

### Pipeline Workflow

```mermaid
flowchart LR
    A[Target ID] --> B[Validation]
    B --> C[Library Prep]
    C --> D[Stage 1: HTVS]
    D --> E[Stage 2: Refined]
    E --> F[Stage 3: Precision]
    F --> G[ADMET Filter]
    G --> H[MD Sims]
    
    style A fill:#22c55e
    style B fill:#22c55e
    style C fill:#22c55e
    style D fill:#22c55e
    style E fill:#f97316
```

---

## 📊 Current Progress

| Phase | Description | Status |
|-------|-------------|--------|
| Phase 0 | Allosteric Target Identification | ✅ Complete |
| Phase 1 | Protocol Validation (AM251 Control) | ✅ Complete |
| Phase 2 | Library Prep & HTVS | 🔄 **Stage 2 Running** |
| Phase 3 | Hit Analysis & Clustering | ⬜ Pending |
| Phase 4 | MD Simulations | ⬜ Pending |
| Phase 5 | Final Analysis & Publication | ⬜ Pending |

### HTVS Pipeline

| Stage | Exhaustiveness | Input | Output | Status |
|-------|----------------|-------|--------|--------|
| Stage 1 | 8 | 10,940 compounds | Top 10% | ✅ Complete |
| Stage 2 | 16 | 2,080 ligands | Top 250/target | 🔄 Running |
| Stage 3 | 32 | 750 ligands | Top 75/target | ⬜ Pending |
| Stage 4 | 64 | 225 ligands | Top 3/target | ⬜ Pending |

---

## 🧬 Target Sites

Three validated binding sites on GPR55:

| Site | Location | Purpose |
|------|----------|---------|
| P0 | Orthosteric | Classical binding pocket |
| P3 | Allosteric | Novel allosteric site |
| Interface | Gα PPI | Protein-protein interface |

---

## 🏆 Top Hits (Stage 2 - In Progress)

| Compound | Target | Affinity |
|----------|--------|----------|
| compound_11569386_3d | P0 | -12.3 kcal/mol |
| compound_121036492_3d | P0 | -11.5 kcal/mol |
| compound_149348105_3d | P0 | -11.4 kcal/mol |

---

## 🛠️ Technologies

- **Docking:** AutoDock Vina 1.2.5
- **Visualization:** PyMOL, UCSF Chimera
- **Cheminformatics:** RDKit, Open Babel
- **Compute:** Google Cloud VM (8-core)

---

## 📁 Repository Structure

```
├── docs/                    # Project documentation
│   ├── GPR55_INTRODUCTION.md
│   ├── PHASE_0_1_COMPLETED.md
│   ├── PHASE_2_PREPARATION.md
│   └── TARGET_COORDINATES.md
├── scripts/
│   ├── colab/              # Colab notebook scripts
│   ├── converters/         # SMILES/SDF/PDBQT converters
│   └── vm/                 # VM docking scripts
├── config/                 # Docking configuration files
├── results/                # Control validation results
└── figures/                # Visualizations & renders
```

---

## 📖 Key Documentation

- [GPR55 Introduction](docs/GPR55_INTRODUCTION.md) - Background on the target
- [Phase 0-1 Complete](docs/PHASE_0_1_COMPLETED.md) - Validation results
- [Phase 2 Preparation](docs/PHASE_2_PREPARATION.md) - Current HTVS setup
- [Target Coordinates](docs/TARGET_COORDINATES.md) - Binding site definitions

---

## 📚 References

1. Eberhardt, J., et al. (2021). AutoDock Vina 1.2.0. *J. Chem. Inf. Model.* DOI: 10.1021/acs.jcim.1c00203
2. Lauckner, J.E., et al. (2008). GPR55 is a cannabinoid receptor. *PNAS* 105(7):2699-704

---

## 📜 License

This project is licensed under the MIT License - see [LICENSE](LICENSE) for details.

---

## 👤 Author

**Ayush** - [GitHub](https://github.com/Aayush-ob)

---

*Part of Bioinformatics Major Project | Timeline: Nov 2025 - Jun 2026*
