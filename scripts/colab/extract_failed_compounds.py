"""
Failed & Unprocessed Compounds Extractor
========================================
Compares original TSV files with conversion results to identify:
- Explicitly failed molecules (from FAILED.tsv)
- Unprocessed molecules (interrupted processing)

Creates a combined file for retry processing.

Run: python extract_failed_compounds.py
"""

import pandas as pd
from rdkit import Chem
import os
from datetime import datetime


# ============================================================================
# CONFIGURATION
# ============================================================================
INPUT_FILES = [
    "GPCR_ligand.tsv",
    "human_protein_ligand.tsv"
]

SMILES_COLUMN = "Smiles"
OUTPUT_FILE = "combined_failed_and_remaining.tsv"
REPORT_FILE = "processing_status_report.txt"


def format_number(num):
    """Format large numbers with commas"""
    return f"{num:,}"


def extract_successful_ids_from_sdf(sdf_file):
    """Extract molecule IDs from an SDF file"""
    successful_ids = set()

    if not os.path.exists(sdf_file):
        print(f"   ⚠️  SDF file not found: {sdf_file}")
        return successful_ids

    try:
        suppl = Chem.SDMolSupplier(sdf_file)
        for mol in suppl:
            if mol is not None:
                mol_id = mol.GetProp("_Name") if mol.HasProp("_Name") else None
                if mol_id:
                    successful_ids.add(mol_id)
        print(f"   ✓ Found {format_number(len(successful_ids))} successful molecules in {sdf_file}")
    except Exception as e:
        print(f"   ✗ Error reading {sdf_file}: {e}")

    return successful_ids


def extract_failed_ids_from_tsv(failed_tsv):
    """Extract molecule IDs from failed TSV file"""
    failed_ids = set()

    if not os.path.exists(failed_tsv):
        print(f"   ⚠️  Failed TSV not found: {failed_tsv}")
        return failed_ids

    try:
        df = pd.read_csv(failed_tsv, sep='\t')
        if 'Molecule_ID' in df.columns:
            failed_ids = set(df['Molecule_ID'].astype(str))
            print(f"   ✓ Found {format_number(len(failed_ids))} explicitly failed molecules in {failed_tsv}")
        else:
            print(f"   ⚠️  No 'Molecule_ID' column in {failed_tsv}")
    except Exception as e:
        print(f"   ✗ Error reading {failed_tsv}: {e}")

    return failed_ids


def get_molecule_id(row, df_columns, idx):
    """Extract molecule ID from a row"""
    for id_col in ['ChEMBL ID', 'Compound ID', 'ID', 'Name', 'Molecule ID']:
        if id_col in df_columns and pd.notna(row.get(id_col)):
            return str(row[id_col])
    return f"Row_{idx}"


def process_file(input_file):
    """Process a single input file and identify failed/unprocessed compounds"""
    base_name = input_file.replace(".tsv", "")
    sdf_file = f"{base_name}_3D_optimized.sdf"
    failed_tsv = f"{base_name}_FAILED.tsv"

    print(f"\n{'='*80}")
    print(f"📂 PROCESSING: {input_file}")
    print('=' * 80)

    # Load original file
    try:
        df_original = pd.read_csv(input_file, sep='\t', low_memory=False)
        print(f"✓ Loaded {format_number(len(df_original))} total rows from {input_file}")
    except Exception as e:
        print(f"✗ ERROR: Could not read {input_file}: {e}")
        return None

    # Check SMILES column
    if SMILES_COLUMN not in df_original.columns:
        print(f"✗ ERROR: Column '{SMILES_COLUMN}' not found in {input_file}")
        return None

    # Remove rows with missing SMILES
    df_original = df_original.dropna(subset=[SMILES_COLUMN])
    total_with_smiles = len(df_original)
    print(f"✓ Found {format_number(total_with_smiles)} rows with valid SMILES")

    # Create molecule ID column
    df_original['Molecule_ID'] = df_original.apply(
        lambda row: get_molecule_id(row, df_original.columns, row.name),
        axis=1
    )

    # Extract successful IDs from SDF
    print(f"\n📊 Checking successful conversions...")
    successful_ids = extract_successful_ids_from_sdf(sdf_file)

    # Extract failed IDs
    print(f"\n📊 Checking explicitly failed molecules...")
    failed_ids = extract_failed_ids_from_tsv(failed_tsv)

    # Identify all processed IDs
    processed_ids = successful_ids.union(failed_ids)

    # Find unprocessed molecules
    all_original_ids = set(df_original['Molecule_ID'])
    unprocessed_ids = all_original_ids - processed_ids

    # Combine failed and unprocessed
    failed_and_unprocessed_ids = failed_ids.union(unprocessed_ids)

    # Filter dataframe
    df_failed_unprocessed = df_original[
        df_original['Molecule_ID'].isin(failed_and_unprocessed_ids)
    ].copy()

    # Add status column
    df_failed_unprocessed['Processing_Status'] = df_failed_unprocessed['Molecule_ID'].apply(
        lambda x: 'Failed_Conversion' if x in failed_ids else 'Not_Processed'
    )

    # Add source file column
    df_failed_unprocessed['Source_File'] = input_file

    # Print summary
    print(f"\n📊 Summary for {input_file}:")
    print(f"   • Total rows with SMILES: {format_number(total_with_smiles)}")
    print(f"   • Successfully converted: {format_number(len(successful_ids))} ({int(len(successful_ids)*100/total_with_smiles) if total_with_smiles > 0 else 0}%)")
    print(f"   • Explicitly failed: {format_number(len(failed_ids))}")
    print(f"   • Not processed (interrupted): {format_number(len(unprocessed_ids))}")
    print(f"   • Total to retry: {format_number(len(failed_and_unprocessed_ids))}")

    return {
        'dataframe': df_failed_unprocessed,
        'stats': {
            'file': input_file,
            'total': total_with_smiles,
            'successful': len(successful_ids),
            'failed': len(failed_ids),
            'unprocessed': len(unprocessed_ids),
            'to_retry': len(failed_and_unprocessed_ids)
        }
    }


def main():
    print("\n" + "=" * 80)
    print(" 🔍 FAILED & UNPROCESSED COMPOUNDS EXTRACTOR")
    print("=" * 80)

    start_time = datetime.now()
    all_results = []
    all_stats = []

    # Process each file
    for input_file in INPUT_FILES:
        if not os.path.exists(input_file):
            print(f"\n⚠️  WARNING: File not found - {input_file}")
            continue

        result = process_file(input_file)
        if result:
            all_results.append(result['dataframe'])
            all_stats.append(result['stats'])

    # Combine results
    if all_results:
        print(f"\n{'='*80}")
        print("📝 COMBINING RESULTS")
        print('=' * 80)

        combined_df = pd.concat(all_results, ignore_index=True)
        combined_df = combined_df.sort_values(['Source_File', 'Processing_Status'])

        # Save combined file
        combined_df.to_csv(OUTPUT_FILE, sep='\t', index=False)
        print(f"✓ Combined file saved: {OUTPUT_FILE}")
        print(f"  Total compounds to retry: {format_number(len(combined_df))}")

        # Generate report
        with open(REPORT_FILE, 'w') as f:
            f.write("=" * 80 + "\n")
            f.write(" PROCESSING STATUS REPORT\n")
            f.write(f" Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
            f.write("=" * 80 + "\n\n")

            total_all = sum(s['total'] for s in all_stats)
            successful_all = sum(s['successful'] for s in all_stats)
            to_retry_all = sum(s['to_retry'] for s in all_stats)

            f.write(f"Total molecules: {format_number(total_all)}\n")
            f.write(f"Successfully converted: {format_number(successful_all)}\n")
            f.write(f"To retry: {format_number(to_retry_all)}\n")

        print(f"✓ Report saved: {REPORT_FILE}")

        # Summary
        elapsed = (datetime.now() - start_time).total_seconds()
        print(f"\n✅ EXTRACTION COMPLETE in {elapsed:.1f}s")
        print(f"\n💡 Next: Run retry converter with {OUTPUT_FILE}")
    else:
        print("\n❌ No results to combine.")


if __name__ == "__main__":
    main()
