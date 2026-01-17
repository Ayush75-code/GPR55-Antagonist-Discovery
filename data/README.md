# Data Directory

This folder contains curated datasets from the GPR55 antagonist virtual screening campaign.

## Files

| File | Description |
|------|-------------|
| `ligand_manifest.csv` | Complete list of all screened ligands with IDs and sources |
| `htvs_results_summary.csv` | Docking results summary (affinities, hit status) |

## For Researchers

### To Download Original Structures
Ligands are sourced from public databases:
- **ChEMBL**: https://www.ebi.ac.uk/chembl/ (use CHEMBL ID)
- **PubChem**: https://pubchem.ncbi.nlm.nih.gov/ (use CID or SID)

### Example: Download CHEMBL5095040
```python
from chembl_webresource_client.new_client import new_client
molecule = new_client.molecule
mol = molecule.get('CHEMBL5095040')
print(mol['molecule_structures']['canonical_smiles'])
```

### Example: Download from PubChem
```python
import pubchempy as pcp
compound = pcp.get_compounds('66844170', 'cid')[0]
print(compound.isomeric_smiles)
```

## Raw Data
Full docking output files (PDBQT, logs) are available upon request.
Contact: [Your Email]
