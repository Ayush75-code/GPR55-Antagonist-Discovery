# Computational Methods & Scripts Documentation

**Project:** GPR55 Antagonist Discovery  
**Date:** 2026-01-17  
**Author:** Ayush Kumar Dewangan

---

## 1. Lipinski Rule of 5 Filtering

### 1.1 Objective
Filter 35,002 compound library to retain only drug-like molecules suitable for oral bioavailability.

### 1.2 Lipinski Criteria
| Parameter | Threshold | Scientific Basis |
|-----------|-----------|------------------|
| Molecular Weight | ≤ 500 Da | Membrane permeability |
| LogP | ≤ 5 | Lipophilicity balance |
| H-bond Donors | ≤ 5 | Intestinal absorption |
| H-bond Acceptors | ≤ 10 | Intestinal absorption |

### 1.3 Script: `filter_lipinski.py`

```python
from rdkit import Chem
from rdkit.Chem import Descriptors, Lipinski
import os, shutil

def check_lipinski(pdbqt_file):
    """Extract SMILES from PDBQT and check Lipinski criteria"""
    with open(pdbqt_file, 'r') as f:
        for line in f:
            if line.startswith("REMARK SMILES"):
                smiles = line.split()[2]
                mol = Chem.MolFromSmiles(smiles)
                if mol:
                    mw = Descriptors.MolWt(mol)
                    logp = Descriptors.MolLogP(mol)
                    hbd = Lipinski.NumHDonors(mol)
                    hba = Lipinski.NumHAcceptors(mol)
                    return mw <= 500 and logp <= 5 and hbd <= 5 and hba <= 10
    return True  # Keep if can't parse
```

### 1.4 Command Line Execution
```powershell
# Install RDKit
pip install rdkit

# Run filtering script
python filter_lipinski.py
```

### 1.5 Results
- **Input:** 35,002 ligands (14 batches)
- **Output:** 11,336 ligands (32% retained)
- **Excluded:** 23,666 compounds (peptides, large molecules)

---

## 2. Interface Target Site Calculation

### 2.1 Objective
Calculate the GPR55-Gα protein interface binding site coordinates for antagonist docking.

### 2.2 Method: PyMOL Interface Analysis

```python
# PyMOL Script: Calculate interface center
# Select interface residues (within 4Å of partner chain)
select interface_residues, (chain R and within 4 of chain A) or (chain A and within 4 of chain R)

# Get center of mass
com = cmd.centerofmass("interface_residues")
print(f"Interface Center: {com}")
# Output: (112.678, 92.605, 112.545)
```

### 2.3 Interface Coordinates
| Parameter | Value |
|-----------|-------|
| Center X | 112.678 Å |
| Center Y | 92.605 Å |
| Center Z | 112.545 Å |
| Box Size | 25 × 25 × 25 Å |

### 2.4 Vina Configuration File: `conf_Interface.txt`
```
center_x = 112.678
center_y = 92.605
center_z = 112.545
size_x = 25
size_y = 25
size_z = 25
```

---

## 3. Batch Reorganization

### 3.1 Script: `reorganize_batches.py`

```python
import os, shutil, math

SOURCE_DIR = "ligand_batches_filtered"
OUTPUT_DIR = "ligands_5batches"

all_files = []
for batch in os.listdir(SOURCE_DIR):
    batch_path = os.path.join(SOURCE_DIR, batch)
    if os.path.isdir(batch_path):
        for f in os.listdir(batch_path):
            if f.endswith('.pdbqt'):
                all_files.append(os.path.join(batch_path, f))

batch_size = math.ceil(len(all_files) / 5)
for i in range(5):
    batch_out = os.path.join(OUTPUT_DIR, f"filtered_batch_{i+1:02d}")
    os.makedirs(batch_out, exist_ok=True)
    for src in all_files[i*batch_size:(i+1)*batch_size]:
        shutil.copy2(src, batch_out)
```

### 3.2 Output Structure
```
ligands_5batches/
├── filtered_batch_01/  (2,268 ligands)
├── filtered_batch_02/  (2,268 ligands)
├── filtered_batch_03/  (2,268 ligands)
├── filtered_batch_04/  (2,268 ligands)
└── filtered_batch_05/  (2,264 ligands)
```

---

## 4. High-Throughput Virtual Screening (HTVS)

### 4.1 AutoDock Vina Parameters
| Parameter | Value | Rationale |
|-----------|-------|-----------|
| Exhaustiveness | 8 | Fast screening mode |
| num_modes | 1 | Quick filter (top pose only) |
| CPU | 1 per job | Parallel efficiency |
| Parallel Jobs | 35-42 | Optimized for 44-core system |

### 4.2 HTVS Script: `run_htvs.sh`

```bash
#!/bin/bash
# AutoDock Vina Parallel HTVS

VINA="/path/to/vina"
RECEPTOR="protein_apo_receptor.pdbqt"
TARGETS=("P0:conf_P0.txt" "P3:conf_P3.txt" "Interface:conf_Interface.txt")

run_vina_job() {
    $VINA --receptor "$RECEPTOR" \
          --ligand "$ligand" \
          --config "$config" \
          --exhaustiveness 8 \
          --num_modes 1 \
          --cpu 1 \
          --out "$output"
}

export -f run_vina_job
cat joblist.txt | parallel --line-buffer -j 35 run_vina_job {}
```

### 4.3 Google Colab Environment
- **Platform:** Google Colab (TPU backend)
- **CPU:** 44-core AMD EPYC 9B14
- **RAM:** 172 GB
- **Storage:** Google Drive mounted

### 4.4 Command Line Execution
```bash
# Install GNU Parallel
apt-get install -y parallel

# Make script executable
chmod +x /tmp/run_htvs.sh

# Run HTVS for batch
/tmp/run_htvs.sh 01
```

---

## 5. Software Versions

| Software | Version | Purpose |
|----------|---------|---------|
| AutoDock Vina | 1.2.5 | Molecular docking |
| RDKit | 2025.9.3 | Cheminformatics |
| GNU Parallel | 20210822 | Parallel processing |
| Python | 3.14 | Scripting |
| PyMOL | 2.5+ | Structure visualization |
| Open Babel | 3.1.1 | File format conversion |

---

## 6. File Locations

| File/Directory | Path | Description |
|----------------|------|-------------|
| Original ligands | `ligand_batches/` | 35,002 compounds (14 batches) |
| Filtered ligands | `ligands_5batches/` | 11,336 compounds (5 batches) |
| Receptor | `Phase 2/protein_apo_receptor.pdbqt` | GPR55 apo structure |
| Config files | `Phase 2/conf_*.txt` | Vina docking parameters |
| Results | `Phase 2/stage1_results/` | Docking output |

---

*Last Updated: 2026-01-17*
