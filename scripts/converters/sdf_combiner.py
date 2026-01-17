#!/usr/bin/env python3
"""
SDF File Combiner - Combines multiple individual SDF files into a single multi-molecule SDF file
"""

import os
import sys
from pathlib import Path


def combine_sdf_files(input_dir, output_file, recursive=False):
    """
    Combine multiple individual SDF files into a single multi-molecule SDF file.
    
    Args:
        input_dir: Directory containing individual SDF files
        output_file: Path to the output combined SDF file
        recursive: If True, search subdirectories for SDF files
    """
    input_path = Path(input_dir)
    output_path = Path(output_file)
    
    if not input_path.exists():
        raise ValueError(f"Input directory '{input_dir}' does not exist!")
    
    if not input_path.is_dir():
        raise ValueError(f"'{input_dir}' is not a directory!")
    
    # Find all SDF files
    if recursive:
        sdf_files = sorted(input_path.rglob('*.sdf'))
    else:
        sdf_files = sorted(input_path.glob('*.sdf'))
    
    if not sdf_files:
        print(f"No SDF files found in {input_dir}")
        return
    
    print(f"Found {len(sdf_files)} SDF files in {input_dir}")
    print(f"Output file: {output_path}")
    print("-" * 60)
    
    # Combine all SDF files
    combined_count = 0
    
    with open(output_path, 'w') as outfile:
        for idx, sdf_file in enumerate(sdf_files, 1):
            try:
                with open(sdf_file, 'r') as infile:
                    content = infile.read().strip()
                    
                    # Skip empty files
                    if not content:
                        print(f"[{idx}/{len(sdf_files)}] Skipped (empty): {sdf_file.name}")
                        continue
                    
                    # Write the content
                    outfile.write(content)
                    
                    # Ensure proper SDF delimiter
                    if not content.endswith('$$$$'):
                        outfile.write('\n$$$$')
                    
                    outfile.write('\n')
                    
                    combined_count += 1
                    
                    if idx % 1000 == 0:
                        print(f"[{idx}/{len(sdf_files)}] Progress: {combined_count} molecules combined")
                    
            except Exception as e:
                print(f"[{idx}/{len(sdf_files)}] Error reading {sdf_file.name}: {e}")
    
    print("-" * 60)
    print(f"✓ Successfully combined {combined_count} molecules")
    print(f"✓ Output file: {output_path.absolute()}")
    print(f"✓ File size: {output_path.stat().st_size / (1024*1024):.2f} MB")


def main():
    """Main function with example usage"""
    if len(sys.argv) < 3:
        print("Usage: python sdf_combiner.py <input_directory> <output_file> [--recursive]")
        print("\nExample:")
        print("  python sdf_combiner.py ./molecules_split combined.sdf")
        print("  python sdf_combiner.py ./molecules_split combined.sdf --recursive")
        print("\nOptions:")
        print("  --recursive, -r  : Search subdirectories for SDF files")
        return
    
    input_dir = sys.argv[1]
    output_file = sys.argv[2]
    recursive = '--recursive' in sys.argv or '-r' in sys.argv
    
    try:
        combine_sdf_files(input_dir, output_file, recursive)
    except Exception as e:
        print(f"Error: {e}")
        import traceback
        traceback.print_exc()


if __name__ == "__main__":
    main()
