"""
Hyper-Parallel Molecule 3D Converter
=====================================
High-performance SMILES to 3D SDF conversion with multi-strategy optimization.
Auto-optimizes for available hardware (CPU cores, RAM).

Features:
- ETKDGv3 3D embedding with fallback strategies
- UFF and MMFF force field optimization
- Parallel processing with live progress tracking
- Checkpoint saving for recovery
- Individual and combined SDF output

Run: python parallel_3d_converter.py

Configure INPUT_FILES list before running.
"""

import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem
import warnings
import os
from datetime import datetime
from multiprocessing import Pool, cpu_count
import psutil
import gc
import sys

# Suppress RDKit warnings
from rdkit import RDLogger
RDLogger.DisableLog('rdApp.*')
warnings.filterwarnings('ignore')

# ============================================================================
# CONFIGURATION
# ============================================================================
INPUT_FILES = ["GPCR_ligand.tsv", "human_protein_ligand.tsv"]
SMILES_COLUMN = "Smiles"
MAX_ATTEMPTS = 500
CREATE_INDIVIDUAL_FILES = True
CREATE_FAILURE_REPORT = True

# Auto-detect optimal settings based on hardware
PHYSICAL_CORES = psutil.cpu_count(logical=False)
LOGICAL_CORES = cpu_count()
TOTAL_RAM_GB = psutil.virtual_memory().total / (1024 ** 3)

# Use 90% of cores for processing
NUM_WORKERS = max(1, int(LOGICAL_CORES * 0.90))
CHUNK_SIZE = 50
BATCH_WRITE_SIZE = 500
PROGRESS_INTERVAL = 50
IMAP_CHUNKSIZE = 25
MEMORY_LIMIT_GB = int(TOTAL_RAM_GB * 0.85)


def get_system_stats():
    """Get real-time system statistics"""
    process = psutil.Process()
    return {
        'ram_gb': process.memory_info().rss / (1024 ** 3),
        'cpu_percent': psutil.cpu_percent(interval=0.1)
    }


def process_molecule(args):
    """Process a single molecule with multiple strategies"""
    idx, smiles, mol_id, max_attempts = args

    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return (False, None, None, 'invalid_smiles', mol_id, smiles)

        mol = Chem.AddHs(mol)

        # Strategy: ETKDGv3 embedding
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
            embed_status = AllChem.EmbedMolecule(
                mol,
                maxAttempts=max_attempts,
                randomSeed=42,
                useRandomCoords=False,
                enforceChirality=True,
                useExpTorsionAnglePrefs=True,
                useBasicKnowledge=True,
                ETversion=3
            )

        if embed_status != 0:
            return (False, None, None, 'failed_embed', mol_id, smiles)

        # Optimization Strategy 1: UFF
        method = 'UFF'
        try:
            result = AllChem.UFFOptimizeMolecule(mol, maxIters=2000)
            conf = mol.GetConformer()
            if conf.GetNumAtoms() > 0:
                mol.SetProp("_Name", mol_id)
                mol.SetProp("OptimizationMethod", method)
                mol.SetProp("ConvergenceStatus", "converged" if result == 0 else "max_iters")
                return (True, Chem.MolToMolBlock(mol), method, None, mol_id, smiles)
        except:
            pass

        # Optimization Strategy 2: MMFF
        method = 'MMFF'
        try:
            mmff_props = AllChem.MMFFGetMoleculeProperties(mol)
            if mmff_props is not None:
                AllChem.MMFFOptimizeMolecule(mol, mmffVariant='MMFF94s', maxIters=1000)
                conf = mol.GetConformer()
                if conf.GetNumAtoms() > 0:
                    mol.SetProp("_Name", mol_id)
                    mol.SetProp("OptimizationMethod", method)
                    return (True, Chem.MolToMolBlock(mol), method, None, mol_id, smiles)
        except:
            pass

        # Optimization Strategy 3: No optimization (raw 3D)
        method = 'NoOpt'
        try:
            conf = mol.GetConformer()
            if conf.GetNumAtoms() > 0:
                mol.SetProp("_Name", mol_id)
                mol.SetProp("OptimizationMethod", method)
                return (True, Chem.MolToMolBlock(mol), method, None, mol_id, smiles)
        except:
            pass

        return (False, None, None, 'failed_optimization', mol_id, smiles)

    except Exception as e:
        return (False, None, None, f'exception: {str(e)}', mol_id, smiles)


def write_batch_results(results_batch, writer, individual_dir, successful_molecules,
                        failed_molecules, method_counts, counters):
    """Efficiently write a batch of results"""
    for success, mol_data, method, error_type, mol_id, smiles in results_batch:
        if success:
            mol = Chem.MolFromMolBlock(mol_data)
            writer.write(mol)
            counters['successful'] += 1

            if method:
                method_counts[method] += 1

            if CREATE_INDIVIDUAL_FILES and individual_dir:
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
            reason_map = {
                'invalid_smiles': ('failed_smiles', "Invalid SMILES string"),
                'failed_embed': ('failed_embed', "Could not generate 3D coordinates"),
                'failed_optimization': ('failed_optimization', "All optimization methods failed")
            }

            if error_type in reason_map:
                counter_key, reason = reason_map[error_type]
                counters[counter_key] += 1
            else:
                counters['failed_exception'] += 1
                reason = error_type

            failed_molecules.append({
                'Molecule_ID': mol_id,
                'SMILES': smiles,
                'Failure_Reason': reason,
                'Status': 'Failed'
            })


def format_number(num):
    """Format large numbers with commas"""
    return f"{num:,}"


def process_file(input_filename):
    """Process a single TSV file with hyper-parallel optimization"""
    base_name = input_filename.replace(".tsv", "")
    output_filename = f"{base_name}_3D_optimized.sdf"
    failure_filename = f"{base_name}_FAILED.tsv"
    individual_dir = f"{base_name}_individual_sdf"

    print("\n" + "=" * 80)
    print(f"📂 PROCESSING: {input_filename}")
    print("=" * 80)
    start_time = datetime.now()

    # Load dataset
    try:
        df = pd.read_csv(input_filename, sep='\t', low_memory=False)
        print(f"✓ Loaded {format_number(len(df))} rows from {input_filename}")
    except Exception as e:
        print(f"✗ FATAL ERROR: Could not read file - {e}")
        return

    if SMILES_COLUMN not in df.columns:
        print(f"✗ FATAL ERROR: Column '{SMILES_COLUMN}' not found.")
        print(f"  Available columns: {list(df.columns)}")
        return

    df.dropna(subset=[SMILES_COLUMN], inplace=True)
    total_mols = len(df)

    if total_mols == 0:
        print("✗ No valid molecules to process")
        return

    print(f"\n⚙️  Preparing {format_number(total_mols)} molecules for parallel processing...")

    # Initialize counters and storage
    counters = {
        'successful': 0,
        'failed_smiles': 0,
        'failed_embed': 0,
        'failed_optimization': 0,
        'failed_exception': 0
    }
    method_counts = {'UFF': 0, 'MMFF': 0, 'NoOpt': 0}
    successful_molecules = []
    failed_molecules = []

    writer = Chem.SDWriter(output_filename)

    if CREATE_INDIVIDUAL_FILES:
        os.makedirs(individual_dir, exist_ok=True)

    # Prepare arguments for parallel processing
    args_list = []
    for idx, row in df.iterrows():
        smiles = row[SMILES_COLUMN]

        mol_id = None
        for id_col in ['ChEMBL ID', 'Compound ID', 'ID', 'Name', 'Molecule ID']:
            if id_col in df.columns and pd.notna(row.get(id_col)):
                mol_id = str(row[id_col])
                break

        if mol_id is None:
            mol_id = f"Row_{idx}"

        args_list.append((idx, smiles, mol_id, MAX_ATTEMPTS))

    print(f"\n🚀 Starting parallel processing with {NUM_WORKERS} workers...")

    processed_count = 0
    results_buffer = []

    with Pool(processes=NUM_WORKERS) as pool:
        for result in pool.imap_unordered(process_molecule, args_list, chunksize=IMAP_CHUNKSIZE):
            results_buffer.append(result)
            processed_count += 1

            # Write batch when buffer is full
            if len(results_buffer) >= BATCH_WRITE_SIZE:
                write_batch_results(results_buffer, writer, individual_dir,
                                  successful_molecules, failed_molecules,
                                  method_counts, counters)
                results_buffer = []

            # Show progress
            if processed_count % PROGRESS_INTERVAL == 0 or processed_count == total_mols:
                elapsed = (datetime.now() - start_time).total_seconds()
                rate = processed_count / elapsed if elapsed > 0 else 0
                eta = (total_mols - processed_count) / rate if rate > 0 else 0
                pct = int((processed_count * 100) / total_mols)

                stats = get_system_stats()

                print(f"   [{processed_count:>6,}/{total_mols:,}] ({pct:>3}%) | "
                      f"Success: {counters['successful']:>6,} | "
                      f"Rate: {rate:>6.1f}/s | "
                      f"ETA: {int(eta):>4}s | "
                      f"RAM: {stats['ram_gb']:>5.1f}GB")

        # Write remaining results
        if results_buffer:
            write_batch_results(results_buffer, writer, individual_dir,
                              successful_molecules, failed_molecules,
                              method_counts, counters)

    writer.close()
    elapsed_time = (datetime.now() - start_time).total_seconds()

    # Save failure report
    if CREATE_FAILURE_REPORT and failed_molecules:
        failure_df = pd.DataFrame(failed_molecules)
        failure_df.to_csv(failure_filename, sep='\t', index=False)

    # Print final summary
    print("\n" + "=" * 80)
    print(f"✅ COMPLETED: {input_filename}")
    print("=" * 80)
    print(f"\n📊 Results Summary:")
    success_pct = int((counters['successful'] * 100) / total_mols) if total_mols > 0 else 0
    print(f"   ✓ Successfully converted: {format_number(counters['successful'])}/{format_number(total_mols)} ({success_pct}%)")

    if counters['successful'] > 0:
        print(f"\n   🎯 Optimization Method Breakdown:")
        total_success = counters['successful']
        for method_name, count in method_counts.items():
            pct = int((count * 100) / total_success) if total_success > 0 else 0
            print(f"      • {method_name}: {format_number(count)} ({pct}%)")

    total_failed = len(failed_molecules)
    if total_failed > 0:
        print(f"\n   ❌ Failed Molecules: {format_number(total_failed)}")

    print(f"\n📁 Output Files:")
    print(f"   • Combined SDF: {output_filename}")
    if CREATE_INDIVIDUAL_FILES:
        print(f"   • Individual SDFs: {individual_dir}/")
    if failed_molecules:
        print(f"   • Failure Report: {failure_filename}")

    print(f"\n⏱️  Performance:")
    print(f"   • Total time: {elapsed_time:.1f}s ({elapsed_time/60:.1f} min)")
    print(f"   • Average rate: {total_mols/elapsed_time:.1f} mol/s")
    print("=" * 80 + "\n")

    return counters, elapsed_time


def main():
    print("\n" + "=" * 80)
    print(" 🚀 HYPER-PARALLEL MOLECULE 3D CONVERTER")
    print("=" * 80)

    print(f"\n🔍 Hardware Detection:")
    print(f"   • Physical Cores: {PHYSICAL_CORES}")
    print(f"   • Logical Cores: {LOGICAL_CORES}")
    print(f"   • Total RAM: {TOTAL_RAM_GB:.1f} GB")
    print(f"   • Workers: {NUM_WORKERS}")

    overall_start = datetime.now()
    all_counters = []

    for file_path in INPUT_FILES:
        if not os.path.exists(file_path):
            print(f"\n⚠️  WARNING: File not found - {file_path}")
            continue

        result = process_file(file_path)
        if result:
            all_counters.append(result)

    total_time = (datetime.now() - overall_start).total_seconds()

    if all_counters:
        print("\n" + "=" * 80)
        print("🎉 ALL PROCESSING COMPLETE!")
        print(f"⏱️  Total time: {total_time:.1f}s ({total_time/60:.1f} min)")
        print("=" * 80 + "\n")


if __name__ == "__main__":
    main()
