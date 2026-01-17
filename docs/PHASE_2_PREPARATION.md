# Phase 2: Library Assembly & High-Throughput Virtual Screening

**Status:** 🔄 IN PROGRESS  
**Duration:** Months 2-3 (Dec 2025 - Feb 2026)

---

## 2.1 Library Curation - ✅ COMPLETED

### Library Composition

| Source | Original Count | After Lipinski Filter | Retained % | Status |
|--------|----------------|----------------------|------------|--------|
| Screening Set (ChEMBL *H. sapiens*) | ~24,006 | 2,203 | 9.2% | ✅ Filtered |
| PubChem Set (Structural Diversity) | ~10,139 | 8,386 | 82.7% | ✅ Filtered |
| GPR55 Focused Set (IC₅₀ < 10 μM) | Included in ChEMBL | - | - | ✅ Filtered |
| AM251 Analogs (Antagonist Control) | ~40 | 37 | 92.5% | ✅ Filtered |
| ML184 Analogs (Agonist Control) | ~315 | 314 | 99.7% | ✅ Filtered |
| Decoy Set (DUD-E Validation) | ~500 | 396 | 79.2% | ✅ Filtered |
| **TOTAL** | **~35,000** | **11,336** | **32.4%** | **Ready** |

> **Note:** ChEMBL Screening Set had low retention (9.2%) due to containing many large peptides and biologics. PubChem and control sets (AM251, ML184) had high retention as they were pre-selected for drug-likeness.

### PDBQT Batching for HTVS

| Batch | Files | Status |
|-------|-------|--------|
| batch_01 - batch_14 | 2,334 each | ✅ Ready |
| batch_15 | 2,326 | ✅ Ready |
| **TOTAL** | **35,002** | **15 batches** |

**Batch Location:** `MAIN DATA/ligand_batches/`  
**Batching Script:** `Scripts/Colab/ligand_batcher.py`

### 2.1.1 Lipinski Rule of 5 Filtering - ✅ COMPLETED (2026-01-17)

**Objective:** Filter library to retain only drug-like compounds suitable for oral bioavailability.

**Lipinski Criteria Applied:**
| Parameter | Threshold | Purpose |
|-----------|-----------|---------|
| Molecular Weight | ≤ 500 Da | Membrane permeability |
| LogP | ≤ 5 | Lipophilicity balance |
| H-bond Donors | ≤ 5 | Absorption |
| H-bond Acceptors | ≤ 10 | Absorption |

**Filtering Results by Original Batch:**

| Batch | Original | Passed Lipinski | Filtered Out |
|-------|----------|-----------------|--------------|
| batch_01 | 2,500 | 370 | 2,130 |
| batch_02 | 2,500 | 327 | 2,173 |
| batch_03 | 2,500 | 335 | 2,165 |
| batch_04 | 2,500 | 155 | 2,345 |
| batch_05 | 2,500 | 213 | 2,287 |
| batch_06 | 2,500 | 121 | 2,379 |
| batch_07 | 2,500 | 122 | 2,378 |
| batch_08 | 2,500 | 87 | 2,413 |
| batch_09 | 2,500 | 216 | 2,284 |
| batch_10 | 2,500 | 1,136 | 1,364 |
| batch_11 | 2,500 | 1,955 | 545 |
| batch_12 | 2,500 | 2,093 | 407 |
| batch_13 | 2,500 | 2,041 | 459 |
| batch_14 | 2,502 | 2,165 | 337 |
| **TOTAL** | **35,002** | **11,336** | **23,666** |

**Filtered Library Reorganization:**

| New Batch | Ligands | Status |
|-----------|---------|--------|
| filtered_batch_01 | 2,268 | ✅ Ready |
| filtered_batch_02 | 2,268 | ✅ Ready |
| filtered_batch_03 | 2,268 | ✅ Ready |
| filtered_batch_04 | 2,268 | ✅ Ready |
| filtered_batch_05 | 2,264 | ✅ Ready |
| **TOTAL** | **11,336** | **5 batches** |

**Output Location:** `ligands_5batches/`  
**Filtering Script:** `filter_lipinski.py` (uses RDKit)

**Notes:**
- Large peptides and non-drug-like molecules (MW > 500 Da, e.g., CHEMBL1076182 with 856 atoms) were excluded
- Batches 10-14 had higher pass rates as they contained more drug-like compounds
- Filtered library reduces docking time by ~68% while retaining all viable drug candidates

### Conversion Pipeline (Google Colab)

```
┌─────────────────────────────────────────────────────────────────────────┐
│                    LIGAND LIBRARY PREPARATION PIPELINE                   │
├─────────────────────────────────────────────────────────────────────────┤
│                                                                          │
│  1. system_info.py          → Check hardware specs & auto-configure      │
│  2. molecule_diagnostics.py → Debug individual molecule issues           │
│  3. parallel_3d_converter.py → Batch SMILES → 3D SDF conversion          │
│  4. extract_failed_compounds.py → Identify failed/unprocessed            │
│  5. retry_converter.py      → Re-process failed compounds                │
│  6. library_builder.py      → Combine all into unified library           │
│                                                                          │
└─────────────────────────────────────────────────────────────────────────┘
```

**Scripts Location:** `Scripts/Colab/`

**Conversion Details:**
- **Method:** RDKit ETKDGv3 embedding + UFF/MMFF optimization
- **Processing:** Parallel (40 cores on Colab Pro+)
- **Rate:** ~100 mol/s on high-end workstations
- **Output:** Combined SDF + individual molecule files

---

## 2.2 Target Preparation - ✅ COMPLETED

### Receptor Files

| Receptor | Chains | File | Status |
|----------|--------|------|--------|
| Apo-Receptor | R | `GPR55_receptor.pdbqt` | ✅ Ready |
| Holo-Complex | R+A+B+C | `GPR55_complex.pdbqt` | ✅ Ready |

### Configuration Files

| Target | Config File | Center (x, y, z) | Status |
|--------|-------------|------------------|--------|
| P0 (Orthosteric) | `conf_target_P0.txt` | (110.02, 111.68, 81.10) | ✅ Ready |
| P1 (Side) | `conf_target_P1.txt` | (122.01, 102.28, 106.69) | ✅ Ready |
| P2 (Interface) | `conf_target_P2.txt` | (97.23, 106.29, 104.49) | ✅ Ready |
| P3 (Allosteric) | `conf_target_P3.txt` | (105.58, 102.05, 107.35) | ✅ Ready |
| P4 (Lower) | `conf_target_P4.txt` | (111.03, 95.74, 101.59) | ✅ Ready |
| P5 (Surface) | `conf_target_P5.txt` | (111.95, 114.75, 103.41) | ✅ Ready |
| **Interface (Holo)** | `conf_target_INTERFACE.txt` | (112.678, 92.605, 112.545) | ✅ Ready |

---

## 2.3 Phase 1 Validation - ✅ COMPLETED

### AM251 Control Docking Results

| Pocket | Best Affinity | Mean | Expected |
|--------|---------------|------|----------|
| **P0 (Orthosteric)** | **-9.458 kcal/mol** | -8.87 | ✅ Strongest |
| P3 (Allosteric) | -7.702 kcal/mol | -7.24 | ✅ Weaker |
| P1-P5 | -6.5 to -7.5 | varies | ✅ Weaker |

**Validation Scripts:** `Scripts/Colab/vina_*.py`

**Key Metrics:**
- Energy Gap (P0 - P3): **~1.76 kcal/mol**
- Hit Threshold: **ΔG < -8.5 kcal/mol**
- Selectivity Filter: **ΔΔG > 1.5 kcal/mol vs P0**

---

## 2.4 HTVS Execution - ⬜ PENDING

### Planned Screening Strategy

**Set A (Apo-Receptor):** P0, P3, P2  
**Set B (Holo-Complex):** Interface, P4, P5

### AutoDock Vina Parameters

| Parameter | Value |
|-----------|-------|
| Software | AutoDock Vina 1.2.5 |
| Exhaustiveness | 16 (HTVS mode) |
| Modes per run | 45 |
| Ligands | ~35,000 |
| Targets | 3 primary (P0, P3, Interface) |
| Total Docking Runs | ~105,000 |

### HTVS Colab Scripts

| Script | Purpose |
|--------|---------|
| `vina_setup.py` | Install Vina, validate files |
| `vina_docking.sh` | Parallel docking (auto-resume) |
| `vina_analysis.py` | Aggregate results, statistics |
| `vina_visualization.py` | Publication plots |

---

## Figures Generated

**Location:** `Figures/PyMOL_Renders/`

### Phase 1 Validation Figures

| Figure | Description |
|--------|-------------|
| `Fig1_Entire_Complex.png` | Full GPR55-G protein complex |
| `Fig2_GPR55_Only.png` | GPR55 receptor with transparent surface |
| `Fig3_Six_Target_Pockets.png` | All 6 binding pockets (P0-P5) |
| `AM251_Control_P0_Orthosteric.png` | AM251 at P0 (-9.458 kcal/mol) |
| `AM251_Control_P3_Allosteric.png` | AM251 at P3 (-7.702 kcal/mol) |
| `HTVS_Entire_Complex.png` | 3 HTVS target sites (P0, P3, Interface) |

### Expected HTVS Outputs (⬜ Pending)

| Figure | Description |
|--------|-------------|
| `binding_affinity_matrix.png` | Heatmap of all ligand-target affinities |
| `plot1_binding_comparison.png` | Box plots per target |
| `plot3_validation_bars.png` | Bar chart comparison |

---

## Next Steps

1. [x] Prepare `GPR55_receptor.pdbqt` ✅
2. [x] Convert ligand libraries to 3D SDF ✅
3. [x] Create Colab processing scripts ✅
4. [ ] Convert 3D SDF → PDBQT for docking
5. [ ] Run HTVS on 3 primary targets (P0, P3, Interface)
6. [ ] Aggregate results and identify hits
7. [ ] Begin Phase 3 hit analysis

---

*Last Updated: 2026-01-15*
