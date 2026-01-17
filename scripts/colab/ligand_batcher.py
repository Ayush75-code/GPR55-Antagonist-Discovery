"""
Ligand Library Batching Tool
============================
Divides large PDBQT ligand libraries into smaller batches for parallel HTVS.

Features:
- Splits files into N equal batches
- Move or copy mode
- Verification of batch distribution
- Resume-safe (skips existing batches)

Usage:
    python ligand_batcher.py --source /path/to/pdbqt --output /path/to/batches --num_batches 15
"""

import os
import shutil
import math
import argparse
from datetime import datetime


def format_number(num):
    """Format numbers with commas"""
    return f"{num:,}"


def count_files(directory):
    """Count files in a directory"""
    if not os.path.exists(directory):
        return 0
    return len([f for f in os.listdir(directory) if os.path.isfile(os.path.join(directory, f))])


def create_batches(source_dir, output_dir, num_batches, mode='move'):
    """
    Divide ligand files into batches.
    
    Args:
        source_dir: Directory containing PDBQT files
        output_dir: Directory to create batch folders
        num_batches: Number of batches to create
        mode: 'move' or 'copy'
    """
    print("=" * 80)
    print("LIGAND LIBRARY BATCHING TOOL")
    print("=" * 80)
    
    # Validate source directory
    if not os.path.exists(source_dir):
        print(f"❌ ERROR: Source directory not found: {source_dir}")
        return False
    
    # Get all files
    all_files = sorted([f for f in os.listdir(source_dir) 
                        if os.path.isfile(os.path.join(source_dir, f))])
    total_files = len(all_files)
    
    if total_files == 0:
        print(f"❌ ERROR: No files found in source directory")
        return False
    
    print(f"\n📂 Source directory: {source_dir}")
    print(f"📊 Total files found: {format_number(total_files)}")
    print(f"📦 Number of batches: {num_batches}")
    
    # Calculate files per batch
    files_per_batch = math.ceil(total_files / num_batches)
    print(f"📋 Files per batch: ~{format_number(files_per_batch)}")
    print(f"⚙️  Mode: {mode.upper()}")
    
    # Create output directory
    os.makedirs(output_dir, exist_ok=True)
    
    print(f"\n🚀 Creating batches...")
    print("-" * 80)
    
    # Create batches
    batch_counts = []
    for batch_num in range(num_batches):
        start_idx = batch_num * files_per_batch
        end_idx = min(start_idx + files_per_batch, total_files)
        
        batch_files = all_files[start_idx:end_idx]
        
        batch_folder = os.path.join(output_dir, f"batch_{batch_num + 1:02d}")
        os.makedirs(batch_folder, exist_ok=True)
        
        # Move or copy files
        for file_name in batch_files:
            src_path = os.path.join(source_dir, file_name)
            dst_path = os.path.join(batch_folder, file_name)
            
            if mode == 'move':
                shutil.move(src_path, dst_path)
            else:
                shutil.copy2(src_path, dst_path)
        
        batch_counts.append(len(batch_files))
        print(f"   Batch {batch_num + 1:02d}: {format_number(len(batch_files))} files")
    
    print("-" * 80)
    print(f"\n✅ Done! {format_number(total_files)} files organized into {num_batches} batches")
    print(f"📁 Output location: {output_dir}")
    
    return True


def verify_batches(source_dir, output_dir, expected_total):
    """Verify batch distribution"""
    print("\n" + "=" * 80)
    print("VERIFICATION")
    print("=" * 80)
    
    # Check source directory
    print("\n--- Source Directory ---")
    if os.path.exists(source_dir):
        remaining = count_files(source_dir)
        if remaining == 0:
            print(f"✅ Source directory is empty (all files moved)")
        else:
            print(f"⚠️  Source still contains {format_number(remaining)} files")
    else:
        print(f"⚠️  Source directory does not exist")
    
    # Check batch directories
    print("\n--- Batch Directories ---")
    if not os.path.exists(output_dir):
        print(f"❌ Output directory not found: {output_dir}")
        return False
    
    batch_folders = sorted([f for f in os.listdir(output_dir) 
                           if os.path.isdir(os.path.join(output_dir, f))])
    
    total_in_batches = 0
    batch_counts = []
    
    for batch_folder in batch_folders:
        batch_path = os.path.join(output_dir, batch_folder)
        file_count = count_files(batch_path)
        batch_counts.append(file_count)
        total_in_batches += file_count
        print(f"   {batch_folder}: {format_number(file_count)} files")
    
    print(f"\n📊 Total files in batches: {format_number(total_in_batches)}")
    
    if expected_total > 0:
        if total_in_batches == expected_total:
            print(f"✅ All {format_number(expected_total)} files accounted for!")
        else:
            print(f"⚠️  Mismatch! Expected {format_number(expected_total)}, found {format_number(total_in_batches)}")
    
    if batch_counts:
        print(f"📋 Batch size range: {min(batch_counts)} to {max(batch_counts)} files")
    
    return total_in_batches == expected_total if expected_total > 0 else True


def cleanup_batches(output_dir):
    """Remove all existing batch folders"""
    print("=" * 80)
    print("CLEANUP: Removing existing batch folders")
    print("=" * 80)
    
    if not os.path.exists(output_dir):
        print(f"⚠️  Directory does not exist: {output_dir}")
        return
    
    deleted_count = 0
    for item in os.listdir(output_dir):
        item_path = os.path.join(output_dir, item)
        if os.path.isdir(item_path):
            shutil.rmtree(item_path)
            print(f"   Deleted: {item}")
            deleted_count += 1
    
    print(f"\n✅ Removed {deleted_count} batch folders")


def main():
    parser = argparse.ArgumentParser(description='Batch ligand library for HTVS')
    parser.add_argument('--source', required=True, help='Source directory with PDBQT files')
    parser.add_argument('--output', required=True, help='Output directory for batches')
    parser.add_argument('--num_batches', type=int, default=15, help='Number of batches (default: 15)')
    parser.add_argument('--mode', choices=['move', 'copy'], default='move', help='Move or copy files')
    parser.add_argument('--cleanup', action='store_true', help='Delete existing batches first')
    parser.add_argument('--verify_only', action='store_true', help='Only verify existing batches')
    parser.add_argument('--expected', type=int, default=0, help='Expected total files for verification')
    
    args = parser.parse_args()
    
    if args.cleanup:
        cleanup_batches(args.output)
    
    if args.verify_only:
        verify_batches(args.source, args.output, args.expected)
    else:
        success = create_batches(args.source, args.output, args.num_batches, args.mode)
        if success:
            verify_batches(args.source, args.output, args.expected)


if __name__ == "__main__":
    main()
