# Session Notes & Daily Log

> Quick notes and observations during project execution

---

## 2026-01-15 (Session 3) - PDBQT Batching for HTVS

### Goal
Convert SDF library to PDBQT and batch for parallel HTVS execution

### Actions Taken

1. **Converted 3D SDF → PDBQT** (Google Colab)
   - Input: 35,002 SDF files from unified library
   - Output: 35,002 PDBQT files ready for docking
   - Location: `MAIN DATA/ligand_library_pdbqt/`

2. **Created 15 Batches for Parallel HTVS**
   - Each batch: ~2,334 PDBQT files
   - Location: `MAIN DATA/ligand_batches/batch_01` to `batch_15`
   - Verified: All 35,002 files accounted for

3. **Created Batching Script**
   - `Scripts/Colab/ligand_batcher.py`
   - Features: move/copy mode, verification, cleanup

### Batch Distribution

| Batches | Files per Batch |
|---------|-----------------|
| batch_01 - batch_14 | 2,334 |
| batch_15 | 2,326 |
| **Total** | **35,002** |

---

## 2026-01-15 (Session 2) - Script Organization & PyMOL Visualization

### Goal
Organize Colab scripts and create PyMOL visualization figures

### Actions Taken

1. **Organized Google Colab Scripts** (`Scripts/Colab/`)
   - **Ligand Processing Pipeline (6 scripts):**
     - `system_info.py` - Hardware detection
     - `molecule_diagnostics.py` - Debug conversion issues
     - `parallel_3d_converter.py` - SMILES → 3D SDF
     - `extract_failed_compounds.py` - Identify failed molecules
     - `retry_converter.py` - Re-process failed compounds
     - `library_builder.py` - Build unified library
   
   - **AutoDock Vina Pipeline (4 scripts):**
     - `vina_setup.py` - Install Vina, validate files
     - `vina_docking.sh` - Parallel docking (auto-resume)
     - `vina_analysis.py` - Aggregate results
     - `vina_visualization.py` - Publication plots

2. **Created PyMOL Visualization Scripts** (`Scripts/PyMOL_Utils/`)
   - `fig_P0_orthosteric_binding.pml` - AM251 at P0 (-9.458 kcal/mol)
   - `fig_P3_allosteric_binding.pml` - AM251 at P3 (-7.702 kcal/mol)
   - `fig_HTVS_target_sites.pml` - 3 HTVS targets visualization

3. **Updated Documentation**
   - Added visualization color scheme to `TARGET_COORDINATES.md`
   - Updated `PHASE_2_PREPARATION.md` with Colab pipeline details
   - Created comprehensive `Scripts/Colab/README.md`

4. **Created Image Processing Utility**
   - `add_white_background.py` - Convert transparent PNGs to white background

### Visualization Color Scheme Established

| Element | Color |
|---------|-------|
| GPR55 (Chain R) | `slate` |
| Gα (Chain A) | `green` |
| Gβ (Chain B) | `lightpink` |
| Gγ (Chain C/G) | `paleyellow` |
| Ligand carbons | `orange` |
| Interacting residues | `green` carbons |

---

## 2026-01-15 (Session 1) - Interface Target Definition

### Goal
Recover forgotten target coordinates for Phase 2 HTVS

### Actions Taken

1. **Analyzed existing configuration files** in `Phase 1 calidation (Control)/`
   - Found all 6 pocket coordinates (P0-P5)
   
2. **Calculated new interface center** using PyMOL:
   ```python
   select interface_gpr55, chain R and (byres chain A around 5)
   # 397 atoms selected
   
   select interface_galpha, chain A and (byres chain R around 5)
   # 335 atoms selected
   
   centerofmass interface_gpr55 or interface_galpha
   # Result: [112.678, 92.605, 112.545]
   ```

3. **Created new config file:** `conf_target_INTERFACE.txt`
   - Center: (112.678, 92.605, 112.545)
   - Box: 25×25×25 Å

4. **Created Progress_Documentation folder** with:
   - README.md
   - PHASE_0_1_COMPLETED.md
   - PHASE_2_PREPARATION.md
   - TARGET_COORDINATES.md
   - TIMELINE_TRACKER.md
   - SESSION_NOTES.md (this file)

### Key Findings
- P3 is on the apo-receptor surface (for blocking recruitment)
- Interface site is the true GPR55/Gα contact zone (for disrupting active complex)
- Both are needed for the "Ensemble Screening" strategy

---

*Add new session notes above*
