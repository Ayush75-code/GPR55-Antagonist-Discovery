# Target Coordinates Reference

> Quick reference for all docking target sites in the GPR55 project.

---

## Set A: Apo-Receptor Targets (Prevention Strategy)

*Receptor file: `GPR55_receptor.pdbqt` (Chain R only)*

| Target ID | Site Name | Center (x, y, z) | Box Size | Config File | Strategy |
|-----------|-----------|------------------|----------|-------------|----------|
| **A1 / P0** | Orthosteric | (110.02, 111.68, 81.10) | 22×22×22 Å | `conf_target_P0.txt` | Control/Reference |
| **A2 / P3** | Allosteric PPI | (105.58, 102.05, 107.35) | 20×20×20 Å | `conf_target_P3.txt` | Block G-protein recruitment |
| **A3 / P2** | Interface Region | (97.23, 106.29, 104.49) | 20×20×20 Å | `conf_target_P2.txt` | Surface wedge |

---

## Set B: Holo-Complex Targets (Disruption Strategy)

*Receptor file: `GPR55_complex.pdbqt` (Chains R + A + B + C)*

| Target ID | Site Name | Center (x, y, z) | Box Size | Config File | Strategy |
|-----------|-----------|------------------|----------|-------------|----------|
| **B1 / Interface** | GPR55/Gα Interface | **(112.678, 92.605, 112.545)** | 25×25×25 Å | `conf_target_INTERFACE.txt` | Disrupt active signaling |
| **B2 / P4** | Lower Pocket | (111.03, 95.74, 101.59) | 20×20×20 Å | `conf_target_P4.txt` | Secondary interface |
| **B3 / P5** | Surface Groove | (111.95, 114.75, 103.41) | 20×20×20 Å | `conf_target_P5.txt` | Interface periphery |

---

## Interface Calculation Details

### Rationale
The GPR55-Gα interface represents a protein-protein interaction (PPI) site critical for signal transduction. Compounds binding here would physically block G-protein coupling, functioning as antagonists regardless of orthosteric pocket occupancy.

### Method: Center of Mass Calculation
**Software:** PyMOL 2.5  
**Structure:** PDB 8ZX5 (GPR55-G protein complex, Cryo-EM)  
**Cutoff Distance:** 5 Å (standard for PPI interface definition)

**Step 1: Identify interface residues on GPR55**
```python
# Select GPR55 residues within 5Å of Gα subunit
select interface_gpr55, chain R and (byres chain A around 5)
# Result: 397 atoms (includes backbone and sidechains)
```

**Step 2: Identify interface residues on Gα**
```python
# Select Gα residues within 5Å of GPR55
select interface_galpha, chain A and (byres chain R around 5)  
# Result: 335 atoms
```

**Step 3: Calculate center of mass**
```python
# Combined interface center of mass
centerofmass interface_gpr55 or interface_galpha
# Result: [112.678, 92.605, 112.545]
```

### Final Interface Target Parameters
| Parameter | Value |
|-----------|-------|
| Center X | 112.678 Å |
| Center Y | 92.605 Å |
| Center Z | 112.545 Å |
| Box Size | 25 × 25 × 25 Å |
| Exhaustiveness (HTVS) | 8 |
| Exhaustiveness (Control) | 64 |

### Validation Result
AM251 control docking at Interface site: **-6.1 to -7.0 kcal/mol** (weak binding), confirming this is a novel target distinct from the orthosteric pocket (-9.5 kcal/mol).

---

## Key Interface Residues

### GPR55 (Chain R) - Interface Residues
```
R45, R48, R56, R59, R60, R63, R118, R119, R122, R123, 
R133, R208, R224, R227, R228, R287, R288, R289, R291, 
R292, R293, and surrounding residues...
```

### Gα Subunit (Chain A) - Interface Residues  
```
Key helix residues including LYS225 and surrounding contacts
```

---

## Threshold Criteria

| Metric | Value | Source |
|--------|-------|--------|
| **Hit Threshold** | ΔG < -8.5 kcal/mol | Phase 1 (P3 baseline) |
| **Selectivity Filter** | ΔΔG > 1.5 kcal/mol vs P0 | AM251 energy gap |

---

## Visualization Color Scheme (PyMOL)

> Consistent color scheme used across all figures.

### Protein Chains

| Chain | Component | PyMOL Color | Hex Code |
|-------|-----------|-------------|----------|
| R | GPR55 Receptor | `slate` | #7A90C2 |
| A | Gα-Subunit | `green` | #00FF00 |
| B | Gβ-Subunit | `lightpink` | #FFB6C1 |
| C/G | Gγ-Subunit | `paleyellow` | #FFFFCC |

### Target Site Spheres (HTVS)

| Target | Site Name | PyMOL Color | Sphere Radius |
|--------|-----------|-------------|---------------|
| P0 | Orthosteric | `cyan` | 11 Å |
| P3 | Allosteric PPI | `tv_green` | 10 Å |
| Interface | GPR55/Gα | `hotpink` | 12.5 Å |

### Six Binding Pockets (Phase 1 Control)

| Pocket | Site Name | PyMOL Color |
|--------|-----------|-------------|
| P0 | Orthosteric | `tv_orange` |
| P1 | Side Pocket | `tv_yellow` |
| P2 | Interface Region | `violet` |
| P3 | Allosteric PPI | `tv_green` |
| P4 | Lower Pocket | `cyan` |
| P5 | Surface Groove | `hotpink` |

### Ligand Visualization

| Element | PyMOL Color |
|---------|-------------|
| Carbon | `orange` |
| Oxygen | `red` |
| Nitrogen | `blue` |
| Sulfur | `yellow` |
| Chlorine | `green` |

### Interacting Residues

| Element | PyMOL Color |
|---------|-------------|
| Carbon | `green` |
| Oxygen | `red` |
| Nitrogen | `blue` |
| Sulfur | `yellow` |

---

*Reference Date: 2026-01-15*
