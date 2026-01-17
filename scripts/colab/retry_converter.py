"""
Fast Retry Converter for Failed/Unprocessed Compounds
=====================================================
Processes compounds that failed in the initial conversion with:
- More aggressive embedding strategies
- Extended optimization iterations
- Live progress tracking
- Checkpoint saving

Run: python retry_converter.py

Expects input: combined_failed_and_remaining.tsv (from extract_failed_compounds.py)
"""

import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem
import warnings
import os
from datetime import datetime
from multiprocessing import Pool, cpu_count
import psutil
import sys
import time

# Suppress RDKit warnings
from rdkit import RDLogger
RDLogger.DisableLog('rdApp.*')
warnings.filterwarnings('ignore')

# ============================================================================
# CONFIGURATION
# ============================================================================
INPUT_FILE = "combined_failed_and_remaining.tsv"
SMILES_COLUMN = "Smiles"
MAX_ATTEMPTS = 2000  # Higher for retry
CREATE_INDIVIDUAL_FILES = True
SAVE_INTERVAL = 100

# Auto-detect optimal settings
PHYSICAL_CORES = psutil.cpu_count(logical=False)
LOGICAL_CORES = cpu_count()
TOTAL_RAM_GB = psutil.virtual_memory().total / (1024 ** 3)
NUM_WORKERS = max(1, int(LOGICAL_CORES * 0.90))
IMAP_CHUNKSIZE = 10


def get_system_stats():
    """Get real-time system statistics"""
    process = psutil.Process()
    return {
        'ram_gb': process.memory_info().rss / (1024 ** 3),
        'cpu_percent': psutil.cpu_percent(interval=0)
    }


def format_number(num):
    return f"{num:,}"


def format_time(seconds):
    if seconds < 60:
        return f"{int(seconds)}s"
    elif seconds < 3600:
        return f"{int(seconds/60)}m {int(seconds%60)}s"
    else:
        return f"{int(seconds/3600)}h {int((seconds%3600)/60)}m"


def process_molecule(args):
    """Process molecule with multiple aggressive strategies"""
    idx, smiles, mol_id, max_attempts, total_count = args
    start_time = time.time()

    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            elapsed = time.time() - start_time
            return (False, None, None, 'invalid_smiles', mol_id, smiles, idx, total_count, elapsed)

        mol = Chem.AddHs(mol)

        # Strategy 1: Standard ETKDGv3
        try:
            params = AllChem.ETKDGv3()
            params.randomSeed = 42
            params.enforceChirality = True
            params.useExpTorsionAnglePrefs = True
            params.useBasicKnowledge = True
            params.numThreads = 1
            params.maxIterations = max_attempts
            embed_status = AllChem.EmbedMolecule(mol, params)
        except AttributeError:
            embed_status = -1

        # Strategy 2: Relaxed constraints
        if embed_status != 0:
            mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
            try:
                params = AllChem.ETKDGv3()
                params.randomSeed = 42
                params.enforceChirality = False
                params.useExpTorsionAnglePrefs = False
                params.useBasicKnowledge = True
                params.numThreads = 1
                params.maxIterations = max_attempts
                embed_status = AllChem.EmbedMolecule(mol, params)
            except:
                pass

        # Strategy 3: Random coordinates
        if embed_status != 0:
            mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
            try:
                embed_status = AllChem.EmbedMolecule(
                    mol,
                    maxAttempts=max_attempts,
                    randomSeed=42,
                    useRandomCoords=True
                )
            except:
                pass

        if embed_status != 0:
            elapsed = time.time() - start_time
            return (False, None, None, 'failed_embed', mol_id, smiles, idx, total_count, elapsed)

        # Optimization: UFF (extended iterations)
        method = 'UFF'
        strategy = 'UFF_3000'
        try:
            result = AllChem.UFFOptimizeMolecule(mol, maxIters=3000)
            conf = mol.GetConformer()
            if conf.GetNumAtoms() > 0:
                mol.SetProp("_Name", mol_id)
                mol.SetProp("OptimizationMethod", method)
                mol.SetProp("Strategy", strategy)
                elapsed = time.time() - start_time
                return (True, Chem.MolToMolBlock(mol), method, None, mol_id, smiles, idx, total_count, elapsed)
        except:
            pass

        # MMFF fallback
        method = 'MMFF'
        try:
            mmff_props = AllChem.MMFFGetMoleculeProperties(mol)
            if mmff_props is not None:
                AllChem.MMFFOptimizeMolecule(mol, mmffVariant='MMFF94s', maxIters=2000)
                conf = mol.GetConformer()
                if conf.GetNumAtoms() > 0:
                    mol.SetProp("_Name", mol_id)
                    mol.SetProp("OptimizationMethod", method)
                    elapsed = time.time() - start_time
                    return (True, Chem.MolToMolBlock(mol), method, None, mol_id, smiles, idx, total_count, elapsed)
        except:
            pass

        # No optimization (raw 3D)
        method = 'NoOpt'
        try:
            conf = mol.GetConformer()
            if conf.GetNumAtoms() > 0:
                mol.SetProp("_Name", mol_id)
                mol.SetProp("OptimizationMethod", method)
                elapsed = time.time() - start_time
                return (True, Chem.MolToMolBlock(mol), method, None, mol_id, smiles, idx, total_count, elapsed)
        except:
            pass

        elapsed = time.time() - start_time
        return (False, None, None, 'failed_optimization', mol_id, smiles, idx, total_count, elapsed)

    except Exception as e:
        elapsed = time.time() - start_time
        return (False, None, None, f'exception: {str(e)}', mol_id, smiles, idx, total_count, elapsed)


def main():
    output_filename = "retry_3D_optimized.sdf"
    failure_filename = "retry_FAILED.tsv"
    individual_dir = "retry_individual_sdf"
    success_log = "retry_SUCCESS.tsv"

    print("\n" + "=" * 80)
    print(" 🚀 RETRY CONVERTER FOR FAILED COMPOUNDS")
    print("=" * 80)

    if not os.path.exists(INPUT_FILE):
        print(f"\n❌ ERROR: Input file not found - {INPUT_FILE}")
        print(f"Please run extract_failed_compounds.py first.")
        return

    # Load dataset
    df = pd.read_csv(INPUT_FILE, sep='\t', low_memory=False)
    print(f"✓ Loaded {format_number(len(df))} rows from {INPUT_FILE}")

    if SMILES_COLUMN not in df.columns:
        print(f"✗ Column '{SMILES_COLUMN}' not found.")
        return

    df.dropna(subset=[SMILES_COLUMN], inplace=True)
    total_mols = len(df)

    if total_mols == 0:
        print("✗ No valid molecules to process")
        return

    # Initialize
    counters = {
        'successful': 0,
        'failed_smiles': 0,
        'failed_embed': 0,
        'failed_optimization': 0,
        'failed_exception': 0
    }
    successful_molecules = []
    failed_molecules = []
    processing_times = []

    writer = Chem.SDWriter(output_filename)

    if CREATE_INDIVIDUAL_FILES:
        os.makedirs(individual_dir, exist_ok=True)

    # Prepare arguments
    args_list = []
    for idx, row in df.iterrows():
        smiles = row[SMILES_COLUMN]
        mol_id = str(row.get('Molecule_ID', f"Row_{idx}"))
        args_list.append((idx, smiles, mol_id, MAX_ATTEMPTS, total_mols))

    print(f"\n🚀 Starting parallel retry with {NUM_WORKERS} workers...")
    start_time = datetime.now()

    processed_count = 0
    checkpoint_counter = 0

    with Pool(processes=NUM_WORKERS) as pool:
        for result in pool.imap_unordered(process_molecule, args_list, chunksize=IMAP_CHUNKSIZE):
            success, mol_data, method, error_type, mol_id, smiles, mol_idx, total, proc_time = result

            processed_count += 1
            checkpoint_counter += 1
            processing_times.append(proc_time)

            if success:
                mol = Chem.MolFromMolBlock(mol_data)
                writer.write(mol)
                counters['successful'] += 1

                if CREATE_INDIVIDUAL_FILES:
                    safe_id = "".join(c for c in mol_id if c.isalnum() or c in ('_', '-'))
                    individual_file = os.path.join(individual_dir, f"{safe_id}.sdf")
                    individual_writer = Chem.SDWriter(individual_file)
                    individual_writer.write(mol)
                    individual_writer.close()

                successful_molecules.append({
                    'Molecule_ID': mol_id,
                    'SMILES': smiles,
                    'Optimization_Method': method,
                    'Status': 'Success'
                })
            else:
                if error_type == 'invalid_smiles':
                    counters['failed_smiles'] += 1
                elif error_type == 'failed_embed':
                    counters['failed_embed'] += 1
                elif error_type == 'failed_optimization':
                    counters['failed_optimization'] += 1
                else:
                    counters['failed_exception'] += 1

                failed_molecules.append({
                    'Molecule_ID': mol_id,
                    'SMILES': smiles,
                    'Failure_Reason': error_type,
                    'Status': 'Failed'
                })

            # Progress update
            if processed_count % 50 == 0 or processed_count == total_mols:
                elapsed = (datetime.now() - start_time).total_seconds()
                rate = processed_count / elapsed if elapsed > 0 else 0
                eta = (total_mols - processed_count) / rate if rate > 0 else 0
                pct = (processed_count * 100) / total_mols

                print(f"   [{processed_count:>6,}/{total_mols:,}] ({pct:>5.1f}%) | "
                      f"Success: {counters['successful']:>6,} | "
                      f"Rate: {rate:>6.1f}/s | "
                      f"ETA: {format_time(eta)}", end='\r')

            # Checkpoint
            if checkpoint_counter >= SAVE_INTERVAL:
                writer.flush()
                if successful_molecules:
                    pd.DataFrame(successful_molecules).to_csv(success_log, sep='\t', index=False)
                if failed_molecules:
                    pd.DataFrame(failed_molecules).to_csv(failure_filename, sep='\t', index=False)
                checkpoint_counter = 0

    print()
    writer.close()
    elapsed_time = (datetime.now() - start_time).total_seconds()

    # Save final results
    if successful_molecules:
        pd.DataFrame(successful_molecules).to_csv(success_log, sep='\t', index=False)
    if failed_molecules:
        pd.DataFrame(failed_molecules).to_csv(failure_filename, sep='\t', index=False)

    # Summary
    print("\n" + "=" * 80)
    print("✅ RETRY COMPLETE!")
    print("=" * 80)
    success_pct = (counters['successful'] * 100) / total_mols if total_mols > 0 else 0
    print(f"\n📊 Results:")
    print(f"   ✓ Successfully converted: {format_number(counters['successful'])}/{format_number(total_mols)} ({success_pct:.1f}%)")
    print(f"\n📁 Output Files:")
    print(f"   • Combined SDF: {output_filename}")
    if CREATE_INDIVIDUAL_FILES:
        print(f"   • Individual SDFs: {individual_dir}/")
    if failed_molecules:
        print(f"   • Failure Report: {failure_filename}")
    print(f"\n⏱️  Total time: {format_time(elapsed_time)}")
    print(f"   Average rate: {total_mols/elapsed_time:.1f} mol/s")
    print("=" * 80 + "\n")


if __name__ == "__main__":
    main()
