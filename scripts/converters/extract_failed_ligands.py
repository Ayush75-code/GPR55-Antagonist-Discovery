#!/usr/bin/env python3
"""
Extract failed SDF files from error log and copy to a separate folder
For re-processing with advanced converter
"""

import os
import shutil
from pathlib import Path
import argparse
import re


def parse_error_log(log_file):
    """Extract failed filenames from error log with improved parsing"""
    failed_files = []
    
    with open(log_file, 'r', encoding='utf-8', errors='ignore') as f:
        for line in f:
            # Look for lines containing error indicators or .sdf files
            # Pattern: any line with .sdf followed by | Error:
            if '.sdf' in line and ('Error:' in line or 'error' in line.lower()):
                # Try to extract .sdf filename
                match = re.search(r'([\w\-]+\.sdf)', line)
                if match:
                    filename = match.group(1)
                    failed_files.append(filename)
                    print(f"Found failed file: {filename}")
    
    # Remove duplicates while preserving order
    seen = set()
    unique_failed = []
    for f in failed_files:
        if f not in seen:
            seen.add(f)
            unique_failed.append(f)
    
    return unique_failed


def copy_failed_files(source_dir, failed_files, output_dir):
    """Copy failed files to output directory"""
    
    # Create output directory
    Path(output_dir).mkdir(parents=True, exist_ok=True)
    
    copied = 0
    not_found = 0
    
    print(f"\nCopying files from {source_dir} to {output_dir}...\n")
    
    for filename in failed_files:
        source_path = Path(source_dir) / filename
        dest_path = Path(output_dir) / filename
        
        if source_path.exists():
            shutil.copy2(source_path, dest_path)
            print(f"✓ Copied: {filename}")
            copied += 1
        else:
            print(f"✗ Not found: {filename}")
            not_found += 1
    
    print(f"\n{'='*70}")
    print(f"SUMMARY:")
    print(f"  Total failed files in log: {len(failed_files)}")
    print(f"  Successfully copied: {copied}")
    print(f"  Not found: {not_found}")
    print(f"  Output directory: {output_dir}")
    print(f"{'='*70}")
    
    return copied, not_found


def main():
    parser = argparse.ArgumentParser(
        description='Extract failed SDF files from error log'
    )
    parser.add_argument(
        '-l', '--log',
        required=True,
        help='Error log file path'
    )
    parser.add_argument(
        '-s', '--source',
        required=True,
        help='Source directory containing original SDF files'
    )
    parser.add_argument(
        '-o', '--output',
        default='./ligand_library_sdf_failed',
        help='Output directory for failed files (default: ./ligand_library_sdf_failed)'
    )
    
    args = parser.parse_args()
    
    # Check if log file exists
    if not os.path.isfile(args.log):
        print(f"Error: Log file '{args.log}' does not exist")
        return 1
    
    # Check if source directory exists
    if not os.path.isdir(args.source):
        print(f"Error: Source directory '{args.source}' does not exist")
        return 1
    
    print(f"Parsing error log: {args.log}")
    failed_files = parse_error_log(args.log)
    print(f"\nTotal unique failed files found: {len(failed_files)}")
    
    if len(failed_files) == 0:
        print("\nNo failed files found in log. Please check if:")
        print("  1. The log file contains error entries")
        print("  2. The log file path is correct")
        return 1
    
    # Copy files
    copied, not_found = copy_failed_files(args.source, failed_files, args.output)
    
    if copied > 0:
        print(f"\n✓ Ready to process {copied} failed files with advanced converter")
        print(f"\nNext step:")
        print(f"python sdf_to_pdbqt_advanced.py -i {args.output} -o ./ligand_library_pdbqt_recovered")
    
    return 0


if __name__ == "__main__":
    exit(main())
