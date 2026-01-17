"""
Molecule Diagnosis Tool
=======================
Diagnoses 3D conversion issues for individual molecules and sample files.
Useful for debugging SMILES parsing, embedding, and optimization failures.

Run: python molecule_diagnostics.py
"""

import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import rdBase
import warnings
import os

# Suppress warnings
from rdkit import RDLogger
RDLogger.DisableLog('rdApp.*')
warnings.filterwarnings('ignore')


def diagnose_single_molecule(smiles, mol_id="test_mol"):
    """Detailed diagnosis of a single molecule conversion"""
    print(f"\n{'='*60}")
    print(f"Diagnosing: {mol_id}")
    print(f"SMILES: {smiles}")
    print(f"{'='*60}")

    # Step 1: Parse SMILES
    print("\n[1] Parsing SMILES...")
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        print("❌ FAILED: Invalid SMILES string")
        return False
    print(f"✓ Success: {mol.GetNumAtoms()} atoms, {mol.GetNumBonds()} bonds")

    # Step 2: Add hydrogens
    print("\n[2] Adding hydrogens...")
    mol_h = Chem.AddHs(mol)
    print(f"✓ Success: {mol_h.GetNumAtoms()} atoms (with H)")

    # Step 3: Try embedding
    print("\n[3] Attempting 3D embedding...")

    # Try ETKDGv3
    try:
        params = AllChem.ETKDGv3()
        params.randomSeed = 42
        embed_status = AllChem.EmbedMolecule(mol_h, params)
        print(f"   ETKDGv3 result: {embed_status} ({'SUCCESS' if embed_status == 0 else 'FAILED'})")

        if embed_status == 0:
            # Try optimization
            print("\n[4] Optimizing with UFF...")
            try:
                opt_result = AllChem.UFFOptimizeMolecule(mol_h, maxIters=200)
                print(f"✓ UFF optimization: {opt_result} ({'converged' if opt_result == 0 else 'max iters'})")

                # Check if we have valid coordinates
                conf = mol_h.GetConformer()
                if conf.GetNumAtoms() > 0:
                    print(f"✓ Valid conformer with {conf.GetNumAtoms()} atoms")

                    # Try to write MolBlock
                    try:
                        mol_block = Chem.MolToMolBlock(mol_h)
                        print(f"✓ MolBlock generated ({len(mol_block)} chars)")
                        print("\n🎉 SUCCESS: Full conversion pipeline works!")
                        return True
                    except Exception as e:
                        print(f"❌ Failed to generate MolBlock: {e}")
                        return False
            except Exception as e:
                print(f"❌ UFF optimization failed: {e}")
                return False
        else:
            print("❌ Embedding failed")
            return False

    except AttributeError as e:
        print(f"❌ ETKDGv3 not available (old RDKit?): {e}")

        # Try fallback
        print("\n   Trying fallback embedding method...")
        try:
            embed_status = AllChem.EmbedMolecule(mol_h, randomSeed=42)
            print(f"   Fallback result: {embed_status}")
            if embed_status == 0:
                print("✓ Fallback embedding worked!")
                return True
        except Exception as e:
            print(f"❌ Fallback also failed: {e}")
            return False

    except Exception as e:
        print(f"❌ Unexpected error: {e}")
        return False

    return False


def diagnose_file_sample(filename, smiles_col="Smiles", num_samples=10):
    """Diagnose a sample of molecules from the file"""
    print(f"\n{'#'*60}")
    print(f"DIAGNOSING FILE: {filename}")
    print(f"{'#'*60}")

    try:
        df = pd.read_csv(filename, sep='\t', low_memory=False)
        print(f"✓ Loaded {len(df)} rows")

        if smiles_col not in df.columns:
            print(f"❌ Column '{smiles_col}' not found!")
            print(f"Available columns: {list(df.columns)}")
            return

        # Clean data
        df_clean = df.dropna(subset=[smiles_col])
        print(f"✓ {len(df_clean)} rows with valid SMILES")

        # Sample molecules
        sample_df = df_clean.head(num_samples)

        results = []
        for idx, row in sample_df.iterrows():
            smiles = str(row[smiles_col])

            # Get molecule ID
            mol_id = None
            for id_col in ['ChEMBL ID', 'Compound ID', 'ID', 'Name']:
                if id_col in df.columns and pd.notna(row.get(id_col)):
                    mol_id = str(row[id_col])
                    break
            if mol_id is None:
                mol_id = f"mol_{idx}"

            success = diagnose_single_molecule(smiles, mol_id)
            results.append(success)

        # Summary
        success_count = sum(results)
        print(f"\n{'='*60}")
        print(f"SAMPLE DIAGNOSIS SUMMARY")
        print(f"{'='*60}")
        print(f"Tested: {len(results)} molecules")
        print(f"Success: {success_count} ({success_count*100//len(results)}%)")
        print(f"Failed: {len(results) - success_count}")

        if success_count == 0:
            print("\n⚠️  CRITICAL: 0% success rate suggests:")
            print("   1. RDKit installation issue")
            print("   2. ETKDGv3 not available (update RDKit)")
            print("   3. System/environment problem")
            print("   4. Data format issue")

    except Exception as e:
        print(f"❌ Error reading file: {e}")


def check_rdkit_environment():
    """Check RDKit installation and capabilities"""
    print("\n" + "=" * 60)
    print("RDKit Environment Check")
    print("=" * 60)
    print(f"RDKit version: {rdBase.rdkitVersion}")

    # Test ETKDGv3
    try:
        params = AllChem.ETKDGv3()
        print("✓ ETKDGv3 available")

        # Check available attributes
        available_attrs = [attr for attr in dir(params) if not attr.startswith('_')]
        print(f"✓ Params attributes: {len(available_attrs)} available")

    except AttributeError:
        print("✗ ETKDGv3 NOT available (older RDKit version)")


def main():
    check_rdkit_environment()

    # Quick test with aspirin
    print("\n" + "=" * 60)
    print("QUICK TEST: Simple molecule (aspirin)")
    print("=" * 60)
    test_smiles = "CC(=O)Oc1ccccc1C(=O)O"
    diagnose_single_molecule(test_smiles, "aspirin")

    # Test with files if they exist
    for filename in ["GPCR_ligand.tsv", "human_protein_ligand.tsv"]:
        if os.path.exists(filename):
            diagnose_file_sample(filename, num_samples=5)
        else:
            print(f"\n⚠️  File not found: {filename}")


if __name__ == "__main__":
    main()
