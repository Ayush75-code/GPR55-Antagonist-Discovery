#!/usr/bin/env python3
"""
PDBQT File Combiner - Combines multiple individual PDBQT files into a single multi-molecule PDBQT file
For use with AutoDock-GPU batch processing
"""

import os
import sys
import argparse
from pathlib import Path


def smart_combine_pdbqt_files(input_dir, output_file, recursive=False):
    """
    Smart version: Handles various PDBQT file formats and ensures proper multi-molecule format.
    """
    input_path = Path(input_dir)
    output_path = Path(output_file)
    
    if not input_path.exists():
        raise ValueError(f"Input directory '{input_dir}' does not exist!")
    
    if not input_path.is_dir():
        raise ValueError(f"'{input_dir}' is not a directory!")
    
    # Find all PDBQT files
    if recursive:
        pdbqt_files = sorted(input_path.rglob('*.pdbqt'))
    else:
        pdbqt_files = sorted(input_path.glob('*.pdbqt'))
    
    if not pdbqt_files:
        print(f"No PDBQT files found in {input_dir}")
        return
    
    print(f"Found {len(pdbqt_files)} PDBQT files in {input_dir}")
    print(f"Output file: {output_path}")
    print("-" * 60)
    
    combined_count = 0
    
    with open(output_path, 'w') as outfile:
        for idx, pdbqt_file in enumerate(pdbqt_files, 1):
            try:
                with open(pdbqt_file, 'r') as infile:
                    lines = infile.readlines()
                    
                    # Skip empty files
                    if not lines:
                        print(f"[{idx}/{len(pdbqt_files)}] Skipped (empty): {pdbqt_file.name}")
                        continue
                    
                    # Remove empty lines at start/end
                    while lines and not lines[0].strip():
                        lines.pop(0)
                    while lines and not lines[-1].strip():
                        lines.pop()
                    
                    if not lines:
                        print(f"[{idx}/{len(pdbqt_files)}] Skipped (empty after cleanup): {pdbqt_file.name}")
                        continue
                    
                    # Process the molecule
                    molecule_lines = []
                    has_model_header = lines[0].startswith('MODEL')
                    has_endmdl_footer = lines[-1].startswith('ENDMDL')
                    
                    # If it has MODEL header, skip it (we'll add our own)
                    if has_model_header:
                        molecule_lines = lines[1:]  # Skip MODEL line
                    else:
                        molecule_lines = lines
                    
                    # Remove ENDMDL if present
                    if has_endmdl_footer and molecule_lines:
                        molecule_lines = molecule_lines[:-1]
                    
                    # Remove any trailing empty lines
                    while molecule_lines and not molecule_lines[-1].strip():
                        molecule_lines.pop()
                    
                    if not molecule_lines:
                        print(f"[{idx}/{len(pdbqt_files)}] Skipped (no content after processing): {pdbqt_file.name}")
                        continue
                    
                    # Write to output
                    if combined_count > 0:
                        outfile.write('\n')
                    
                    outfile.write(f'MODEL        {combined_count + 1}\n')
                    outfile.writelines(molecule_lines)
                    outfile.write('ENDMDL\n')
                    
                    combined_count += 1
                    
                    if idx % 1000 == 0:
                        print(f"[{idx}/{len(pdbqt_files)}] Progress: {combined_count} molecules combined")
                    
            except Exception as e:
                print(f"[{idx}/{len(pdbqt_files)}] Error reading {pdbqt_file.name}: {e}")
    
    print("-" * 60)
    print(f"✓ Successfully combined {combined_count} molecules")
    print(f"✓ Output file: {output_path.absolute()}")
    print(f"✓ File size: {output_path.stat().st_size / (1024*1024):.2f} MB")


def main():
    """Main function with example usage"""
    parser = argparse.ArgumentParser(
        description='Combine multiple PDBQT files into a single multi-molecule PDBQT file'
    )
    parser.add_argument('input_dir', help='Directory containing individual PDBQT files')
    parser.add_argument('output_file', help='Path to the output combined PDBQT file')
    parser.add_argument('--recursive', '-r', action='store_true', 
                       help='Search subdirectories for PDBQT files')
    
    if len(sys.argv) < 3:
        parser.print_help()
        print("\nExamples:")
        print("  python pdbqt_combiner.py ./ligands combined.pdbqt")
        print("  python pdbqt_combiner.py ./ligands combined.pdbqt --recursive")
        return
    
    args = parser.parse_args()
    
    try:
        smart_combine_pdbqt_files(args.input_dir, args.output_file, args.recursive)
    except Exception as e:
        print(f"Error: {e}")
        import traceback
        traceback.print_exc()


if __name__ == "__main__":
    main()
