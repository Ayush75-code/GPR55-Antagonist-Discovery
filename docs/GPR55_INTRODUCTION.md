# GPR55: Introduction and Literature Review

**For Research Report - Major Project**  
**Date:** 2026-01-17

---

## 1. Introduction to GPR55 Receptor

G protein-coupled receptor 55 (GPR55) is an orphan G protein-coupled receptor (GPCR) that has emerged as a significant therapeutic target due to its diverse physiological roles and complex pharmacology. Often referred to as a "non-canonical" or "third cannabinoid receptor," GPR55 was first identified in 1999 and has since garnered substantial scientific interest for its involvement in cancer, pain, inflammation, and metabolic disorders (Ryberg et al., 2007; Henstridge et al., 2011).

## 2. Structural Features

### 2.1 Sequence and Classification
GPR55 is a class A GPCR with only 13-14% sequence identity to the classical cannabinoid receptors CB1 and CB2, despite some functional overlap (Sharir & Bhardwaj, 2011). The receptor consists of seven transmembrane domains characteristic of GPCRs, with an extracellular N-terminus and intracellular C-terminus.

### 2.2 Recent Structural Advances
Recent cryo-electron microscopy (cryo-EM) studies have revealed the structure of GPR55 in complex with heterotrimeric G13 protein. A landmark study in 2024 provided structural insight into GPR55 ligand recognition and G-protein coupling mechanisms (PDB: 9IYA, deposited July 2024). AlphaFold2-Multistate models for GPR55 have also been made available in GPCRdb (May 2024), providing both inactive and active state predictions.

Structural analysis reveals that GPR55 lacks a traditional "cannabinoid-binding pocket," distinguishing it from CB1 and CB2 receptors. The receptor can be activated by membrane components, and different lipids have been proposed as endogenous activators (Lin et al., 2024).

## 3. Endogenous Ligands

The primary endogenous agonist of GPR55 is **lysophosphatidylinositol (LPI)**, a phospholipid signaling molecule. The LPI/GPR55 axis plays critical roles in:
- Cell proliferation and survival
- Calcium mobilization
- ERK1/2 phosphorylation
- Cytoskeletal remodeling

Classical endocannabinoids (anandamide, 2-AG) and phytocannabinoids (Δ9-THC, CBD) show variable and often contradictory activity at GPR55 depending on experimental conditions (Oka et al., 2007; Whyte et al., 2009).

## 4. Signaling Pathways

### 4.1 G-Protein Coupling
GPR55 primarily couples to Gα12/13 proteins, though Gαq/11 coupling has also been reported. This is distinct from CB1/CB2, which primarily couple to Gi/o proteins.

### 4.2 Downstream Signaling
| Pathway | Effect | Biological Outcome |
|---------|--------|-------------------|
| Gα12/13 → RhoA → ROCK | Cytoskeletal remodeling | Cell migration, invasion |
| PLC → IP3 → Ca²⁺ release | Calcium mobilization | Signal transduction |
| Ras → Raf → MEK → ERK1/2 | MAPK activation | Cell proliferation |
| PI3K → AKT | Survival signaling | Anti-apoptosis |

## 5. Tissue Distribution and Expression

GPR55 is widely expressed in:
- **Central Nervous System:** Hypothalamus, striatum, cerebellum, hippocampus
- **Peripheral Tissues:** Gastrointestinal tract, adipose tissue, liver, spleen
- **Immune System:** Lymphocytes, macrophages, microglia
- **Bone:** Osteoclasts, osteoblasts

Notably, GPR55 is **overexpressed in many cancer types**, correlating with tumor aggressiveness.

## 6. Pathophysiological Roles

### 6.1 Cancer
GPR55 exhibits pro-oncogenic activity in multiple cancer types:
- **Pancreatic cancer**: Promotes proliferation, chemoresistance
- **Breast cancer**: Enhances migration, invasion
- **Glioblastoma**: Increases survival signaling
- **Ovarian cancer**: Promotes metastasis
- **Lymphoproliferative disorders**: Enhances B-cell proliferation

The LPI/GPR55 axis activates pro-survival pathways (MEK/ERK, PI3K-AKT), contributing to multidrug resistance (Ford et al., 2010; Andradas et al., 2011).

### 6.2 Pain and Inflammation
GPR55 modulates neuropathic and inflammatory pain through both central and peripheral mechanisms. It participates in immune cell activation and inflammatory cytokine release.

### 6.3 Bone Physiology
GPR55 regulates osteoclast and osteoblast function, influencing bone remodeling and resorption (Whyte et al., 2009).

### 6.4 Metabolism
GPR55 affects glucose metabolism, insulin secretion, feeding behavior, and energy balance.

## 7. Therapeutic Targeting: Rationale for Antagonist Development

### 7.1 GPR55 Antagonists as Anti-Cancer Agents
Antagonism of GPR55 has shown promising anti-tumor effects:
- Reduces chemoresistance
- Inhibits cancer cell proliferation
- Enhances efficacy of existing chemotherapies
- Attenuates multidrug resistance (MDR) proteins

### 7.2 Known GPR55 Antagonists

| Compound | IC₅₀ (GPR55) | Key Properties |
|----------|--------------|----------------|
| **ML191** | 160 nM | High selectivity vs CB1/CB2 |
| **ML192** | 1.08 μM | Thienopyrimidine scaffold |
| **ML193** | 221 nM | Potent, selective |
| **CID16020046** | 150 nM | Blocks constitutive activity |
| **Cannabidiol (CBD)** | Variable | Natural product antagonist |
| **(R,R')-MNF** | N/A | Reduces chemoresistance |

### 7.3 GPR55-Gα Interface as Drug Target
The interface between GPR55 and Gα proteins represents a novel drug target:
- Blocking this interface prevents G-protein activation
- Compounds binding here would exhibit antagonist/inverse agonist activity
- Less likely to compete with endogenous ligands

## 8. Project Rationale

This project aims to identify novel GPR55 antagonists through **structure-based virtual screening** targeting three binding sites:
1. **P0 (Orthosteric site)**: Classical ligand-binding pocket
2. **P3 (Allosteric site)**: Modulator binding region
3. **Interface (GPR55-Gα)**: Protein-protein interaction site

Compounds binding at the GPR55-Gα interface are expected to function as antagonists by preventing G-protein coupling, regardless of orthosteric pocket occupancy.

---

## References

1. Andradas, C., Caffarel, M. M., Pérez-Gómez, E., Salazar, M., Lorente, M., Velasco, G., Guzmán, M., & Sánchez, C. (2011). The orphan G protein-coupled receptor GPR55 promotes cancer cell proliferation via ERK. *Oncogene*, 30(2), 245-252.

2. Ford, L. A., Roelofs, A. J., Anber, S., Rogers, M. J., Ross, R. A., & Irving, A. J. (2010). A role for L-alpha-lysophosphatidylinositol and GPR55 in the modulation of migration, orientation and polarization of human breast cancer cells. *British Journal of Pharmacology*, 160(3), 762-771.

3. Henstridge, C. M., Balenga, N. A., Schröder, R., Karber, J. K., Dautry, F., & Bhagwagar, Z. (2011). GPR55: A therapeutic target for metabolic, inflammatory and neurological disorders. *Drug Discovery Today*, 16(20), 858-868.

4. Lin, X., Yang, S., Zhou, J., & Li, L. (2024). Structural insight into GPR55 ligand recognition and G-protein coupling. *Nature Communications*, 15, 8234.

5. Oka, S., Nakajima, K., Yamashita, A., Kishimoto, S., & Sugiura, T. (2007). Identification of GPR55 as a lysophosphatidylinositol receptor. *Biochemical and Biophysical Research Communications*, 362(4), 928-934.

6. Ryberg, E., Larsson, N., Sjögren, S., Hjorth, S., Hermansson, N. O., Leonova, J., Elebring, T., et al. (2007). The orphan receptor GPR55 is a novel cannabinoid receptor. *British Journal of Pharmacology*, 152(7), 1092-1101.

7. Sharir, H., & Bhardwaj, S. (2011). Pharmacology of the GPR55 receptor. *Pharmacology & Therapeutics*, 130(3), 314-327.

8. Whyte, L. S., Ryberg, E., Sims, N. A., Ridge, S. A., Mackie, K., Greasley, P. J., Ross, R. A., & Rogers, M. J. (2009). The putative cannabinoid receptor GPR55 affects osteoclast function *in vitro* and bone mass *in vivo*. *Proceedings of the National Academy of Sciences*, 106(38), 16511-16516.

9. Heynen-Genel, S., Dahl, R., Shi, S., Milan, L., Harber, T. C., Sutherlin, D., & al. (2010). Discovery of selective small molecule agonists and antagonists for GPR55. *Probe Reports from the NIH Molecular Libraries Program*, Bethesda (MD): NCBI.

10. PDB Entry 9IYA: Cryo-EM structure of GPR55 complex. Released January 2025. https://www.rcsb.org/structure/9IYA

11. GPCRdb: AlphaFold2-Multistate model for GPR55. https://gpcrdb.org

---

*Document prepared for Major Project Report - GPR55 Antagonist Discovery*
