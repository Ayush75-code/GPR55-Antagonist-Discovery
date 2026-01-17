"""
Ligand Library Builder
======================
Combines all successfully converted molecules from multiple sources
into a unified library with:
- Combined SDF file (all molecules)
- Individual SDF files (one per molecule)
- Molecule catalog (TSV with properties)
- Statistics and documentation

Run: python library_builder.py
"""

import os
import shutil
from pathlib import Path
from rdkit import Chem
from rdkit.Chem import Descriptors
import pandas as pd
from datetime import datetime
from collections import defaultdict


# ============================================================================
# CONFIGURATION
# ============================================================================
# Define source directories (individual SDF folders)
SOURCES = {
    'human_protein_ligand': 'human_protein_ligand_individual_sdf',
    'GPCR_ligand': 'GPCR_ligand_individual_sdf',
    'retry': 'retry_individual_sdf'
}

LIBRARY_DIR = "LIGAND_LIBRARY"


def format_number(num):
    """Format large numbers with commas"""
    return f"{num:,}"


def check_individual_dir(dir_path, source_name):
    """Check individual SDF directory and load all molecules"""
    if not os.path.exists(dir_path):
        print(f"   ✗ Directory not found: {dir_path}")
        return []

    sdf_files = [f for f in os.listdir(dir_path) if f.endswith('.sdf')]
    print(f"   📁 Found {format_number(len(sdf_files))} SDF files in {dir_path}")

    valid_mols = []
    invalid_count = 0

    for i, sdf_file in enumerate(sdf_files, 1):
        file_path = os.path.join(dir_path, sdf_file)
        try:
            suppl = Chem.SDMolSupplier(file_path)
            mol = suppl[0]
            if mol is not None:
                mol_id = mol.GetProp("_Name") if mol.HasProp("_Name") else os.path.splitext(sdf_file)[0]
                valid_mols.append((mol_id, mol, sdf_file))
            else:
                invalid_count += 1
        except Exception:
            invalid_count += 1

        if i % 1000 == 0:
            print(f"      Progress: {i}/{len(sdf_files)}...", end='\r')

    print(f"      ✓ Valid molecules: {format_number(len(valid_mols))}")
    if invalid_count > 0:
        print(f"      ⚠️  Invalid/corrupt files: {format_number(invalid_count)}")

    return valid_mols


def main():
    print("\n" + "=" * 100)
    print("📚 LIGAND LIBRARY BUILDER")
    print("=" * 100)

    # Step 1: Load molecules
    print("\n📂 Step 1: Loading molecules from individual SDF folders...")
    print("-" * 100)

    all_molecules = {}
    source_counts = defaultdict(int)

    for source_name, dir_path in SOURCES.items():
        print(f"\n🔍 Loading from: {source_name}")
        print(f"   Directory: {dir_path}")

        mols = check_individual_dir(dir_path, source_name)

        for mol_id, mol, filename in mols:
            if mol_id not in all_molecules:
                all_molecules[mol_id] = (mol, source_name, filename)
                source_counts[source_name] += 1

    total_unique = len(all_molecules)
    print(f"\n✨ Total Unique Molecules: {format_number(total_unique)}")

    if total_unique == 0:
        print("\n❌ No molecules found!")
        return

    print(f"\n📈 Breakdown by source:")
    for source, count in sorted(source_counts.items(), key=lambda x: x[1], reverse=True):
        pct = (count * 100) / total_unique
        print(f"   • {source}: {format_number(count)} ({pct:.1f}%)")

    # Step 2: Create library directory
    print(f"\n🏗️  Step 2: Creating library directory...")

    if os.path.exists(LIBRARY_DIR):
        print(f"   ⚠️  Directory '{LIBRARY_DIR}' already exists!")
        response = input("   Overwrite? (yes/no): ").strip().lower()
        if response != 'yes':
            print("   ❌ Cancelled.")
            return
        shutil.rmtree(LIBRARY_DIR)

    os.makedirs(LIBRARY_DIR)
    for subdir in ['combined_sdf', 'individual_molecules', 'metadata']:
        os.makedirs(os.path.join(LIBRARY_DIR, subdir))
    print(f"   ✓ Created: {LIBRARY_DIR}")

    # Step 3: Write combined SDF
    print(f"\n💾 Step 3: Writing combined SDF file...")

    combined_sdf = os.path.join(LIBRARY_DIR, 'combined_sdf', 'ALL_LIGANDS_3D.sdf')
    writer = Chem.SDWriter(combined_sdf)

    molecules_list = []
    for mol_id, (mol, source, orig_file) in sorted(all_molecules.items()):
        mol.SetProp("Source", source)
        mol.SetProp("Original_Filename", orig_file)
        mol.SetProp("_Name", mol_id)
        writer.write(mol)

        molecules_list.append({
            'Molecule_ID': mol_id,
            'Source': source,
            'Original_File': orig_file,
            'Has_3D_Coords': mol.GetNumConformers() > 0,
            'Num_Atoms': mol.GetNumAtoms(),
            'Num_Heavy_Atoms': mol.GetNumHeavyAtoms(),
            'Molecular_Weight': round(Descriptors.MolWt(mol), 2)
        })

    writer.close()
    file_size_mb = os.path.getsize(combined_sdf) / (1024 * 1024)
    print(f"   ✓ Written: {format_number(total_unique)} molecules ({file_size_mb:.2f} MB)")

    # Step 4: Write individual files
    print(f"\n📁 Step 4: Copying individual SDF files...")

    individual_dir = os.path.join(LIBRARY_DIR, 'individual_molecules')
    for mol_id, (mol, source, orig_file) in sorted(all_molecules.items()):
        safe_id = "".join(c if c.isalnum() or c in ('_', '-') else '_' for c in mol_id)
        dest_file = os.path.join(individual_dir, f"{safe_id}.sdf")
        individual_writer = Chem.SDWriter(dest_file)
        individual_writer.write(mol)
        individual_writer.close()
    print(f"   ✓ Copied {format_number(total_unique)} individual files")

    # Step 5: Create catalog
    print(f"\n📝 Step 5: Creating molecule catalog...")

    catalog_file = os.path.join(LIBRARY_DIR, 'metadata', 'molecule_catalog.tsv')
    df_catalog = pd.DataFrame(molecules_list)
    df_catalog = df_catalog.sort_values('Molecule_ID')
    df_catalog.to_csv(catalog_file, sep='\t', index=False)
    print(f"   ✓ Created: {catalog_file}")

    # Step 6: Create statistics
    stats_file = os.path.join(LIBRARY_DIR, 'metadata', 'library_statistics.txt')
    with open(stats_file, 'w') as f:
        f.write("=" * 80 + "\n")
        f.write("LIGAND LIBRARY STATISTICS\n")
        f.write("=" * 80 + "\n\n")
        f.write(f"Creation Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
        f.write(f"Total Unique Molecules: {format_number(total_unique)}\n\n")
        for source, count in sorted(source_counts.items(), key=lambda x: x[1], reverse=True):
            f.write(f"  {source}: {format_number(count)}\n")
    print(f"   ✓ Created: {stats_file}")

    # Step 7: Create README
    readme_file = os.path.join(LIBRARY_DIR, 'README.md')
    with open(readme_file, 'w') as f:
        f.write("# Ligand Library\n\n")
        f.write(f"**Created:** {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")
        f.write(f"**Total Molecules:** {format_number(total_unique)}\n\n")
        f.write("## Sources\n\n")
        for source, count in sorted(source_counts.items(), key=lambda x: x[1], reverse=True):
            f.write(f"- **{source}**: {format_number(count)} molecules\n")
        f.write("\n## Quick Usage\n\n```python\nfrom rdkit import Chem\n")
        f.write("suppl = Chem.SDMolSupplier('combined_sdf/ALL_LIGANDS_3D.sdf')\nfor mol in suppl:\n")
        f.write("    if mol: print(mol.GetProp('_Name'))\n```\n")
    print(f"   ✓ Created: {readme_file}")

    # Final summary
    print("\n" + "=" * 100)
    print("✅ LIGAND LIBRARY CREATED!")
    print("=" * 100)
    print(f"\n📚 Location: ./{LIBRARY_DIR}/")
    print(f"   • Total molecules: {format_number(total_unique)}")
    print(f"   • Combined SDF: {file_size_mb:.2f} MB")
    print("\n" + "=" * 100 + "\n")


if __name__ == "__main__":
    main()
