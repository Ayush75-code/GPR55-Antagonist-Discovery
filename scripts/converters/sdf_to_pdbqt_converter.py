#!/usr/bin/env python3
"""
SDF to PDBQT Converter for AutoDock Vina and AutoDock-GPU
Converts all SDF files in a directory to PDBQT format using RDKit + Meeko
Features: Checkpointing, resume capability, separate logs for success/error
"""

import os
import sys
import argparse
import json
from pathlib import Path
from multiprocessing import Pool, cpu_count
from datetime import datetime

# Create logs directory in script location
SCRIPT_DIR = Path(__file__).parent.absolute()
LOGS_DIR = SCRIPT_DIR / "conversion_logs"
LOGS_DIR.mkdir(exist_ok=True)


class Logger:
    """Logger that writes to both console and multiple log files"""
    def __init__(self, main_log, success_log, error_log):
        self.main_log = main_log
        self.success_log = success_log
        self.error_log = error_log
        
        # Open all log files
        self.main_handle = open(main_log, 'a')
        self.success_handle = open(success_log, 'a')
        self.error_handle = open(error_log, 'a')
        
        # Write headers
        timestamp = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
        header = f"=== SDF to PDBQT Conversion Log ===\nSession started at: {timestamp}\n"
        
        self.main_handle.write(header)
        self.main_handle.flush()
    
    def log(self, message):
        """Write message to console and main log file"""
        print(message)
        self.main_handle.write(message + '\n')
        self.main_handle.flush()
    
    def log_success(self, filename, output_file):
        """Log successful conversion"""
        msg = f"✓ {filename} -> {output_file}"
        self.success_handle.write(f"{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} | {msg}\n")
        self.success_handle.flush()
    
    def log_error(self, filename, error):
        """Log failed conversion"""
        msg = f"✗ {filename} | Error: {error}"
        self.error_handle.write(f"{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} | {msg}\n")
        self.error_handle.flush()
    
    def close(self):
        """Close all log files"""
        footer = f"\nSession finished at: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n{'='*70}\n"
        self.main_handle.write(footer)
        
        self.main_handle.close()
        self.success_handle.close()
        self.error_handle.close()


class CheckpointManager:
    """Manages checkpoint file for resume capability"""
    def __init__(self, checkpoint_file):
        self.checkpoint_file = checkpoint_file
        self.processed_files = set()
        self.load()
    
    def load(self):
        """Load checkpoint data"""
        if self.checkpoint_file.exists():
            try:
                with open(self.checkpoint_file, 'r') as f:
                    data = json.load(f)
                    self.processed_files = set(data.get('processed', []))
            except Exception as e:
                print(f"Warning: Could not load checkpoint: {e}")
    
    def save(self, filename, success):
        """Save processed file to checkpoint"""
        self.processed_files.add(filename)
        try:
            with open(self.checkpoint_file, 'w') as f:
                json.dump({'processed': list(self.processed_files)}, f)
        except Exception as e:
            print(f"Warning: Could not save checkpoint: {e}")
    
    def is_processed(self, filename):
        """Check if file was already processed"""
        return filename in self.processed_files
    
    def get_progress(self, total):
        """Get progress statistics"""
        return len(self.processed_files), total


def convert_single_file_meeko(args):
    """Convert a single SDF file to PDBQT"""
    sdf_file, output_dir = args
    
    try:
        from rdkit import Chem
        from rdkit.Chem import AllChem
        from meeko import MoleculePreparation
        from meeko import PDBQTWriterLegacy
        
        # Read molecule
        supplier = Chem.SDMolSupplier(str(sdf_file), removeHs=False, sanitize=True)
        results = []
        
        for idx, mol in enumerate(supplier):
            if mol is None:
                results.append((sdf_file.name, False, "Failed to read molecule"))
                continue
            
            try:
                # Add explicit hydrogens
                mol = Chem.AddHs(mol, addCoords=True)
                
                # Generate 3D coords if missing
                if mol.GetNumConformers() == 0:
                    result = AllChem.EmbedMolecule(mol, randomSeed=42)
                    if result != 0:
                        results.append((sdf_file.name, False, "Failed to generate 3D coordinates"))
                        continue
                    AllChem.MMFFOptimizeMolecule(mol, maxIters=200)
                
                # Prepare with Meeko
                preparator = MoleculePreparation()
                mol_prepared_list = preparator.prepare(mol)
                
                if not mol_prepared_list or len(mol_prepared_list) == 0:
                    results.append((sdf_file.name, False, "Meeko returned empty result"))
                    continue
                
                # Get first prepared molecule
                mol_prepared = mol_prepared_list[0]
                
                # Determine output filename
                num_mols = len([m for m in Chem.SDMolSupplier(str(sdf_file), removeHs=False) if m is not None])
                if num_mols > 1:
                    output_file = Path(output_dir) / f"{sdf_file.stem}_{idx+1}.pdbqt"
                else:
                    output_file = Path(output_dir) / f"{sdf_file.stem}.pdbqt"
                
                # Write PDBQT - handle both string and tuple returns
                pdbqt_output = PDBQTWriterLegacy.write_string(mol_prepared)
                
                # Handle tuple return (pdbqt_string, is_ok) or just string
                if isinstance(pdbqt_output, tuple):
                    pdbqt_string = pdbqt_output[0]
                else:
                    pdbqt_string = pdbqt_output
                
                if not pdbqt_string or len(pdbqt_string) == 0:
                    results.append((sdf_file.name, False, "Empty PDBQT output"))
                    continue
                
                with open(output_file, 'w') as f:
                    f.write(pdbqt_string)
                
                results.append((output_file.name, True, None))
                
            except Exception as e:
                results.append((sdf_file.name, False, f"Processing error: {str(e)}"))
                continue
        
        if not results:
            results.append((sdf_file.name, False, "No valid molecules in file"))
        
        return results
        
    except Exception as e:
        return [(sdf_file.name, False, f"File error: {str(e)}")]


def convert_with_meeko(input_dir, output_dir, num_cores, logger, checkpoint_mgr):
    """Convert SDF files to PDBQT with checkpointing"""
    try:
        from rdkit import Chem
        from meeko import MoleculePreparation
    except ImportError as e:
        msg = f"Error: Required packages not installed.\nImport error: {e}\nInstall: conda install -c conda-forge rdkit && pip install meeko gemmi"
        logger.log(msg)
        sys.exit(1)
    
    logger.log("Using RDKit + Meeko for conversion with checkpointing enabled...")
    
    # Get all SDF files
    all_sdf_files = sorted(list(Path(input_dir).glob("*.sdf")))
    
    if not all_sdf_files:
        logger.log(f"No SDF files found in {input_dir}")
        return
    
    # Filter out already processed files
    sdf_files = [f for f in all_sdf_files if not checkpoint_mgr.is_processed(f.name)]
    
    already_processed = len(all_sdf_files) - len(sdf_files)
    if already_processed > 0:
        logger.log(f"Resuming: {already_processed} files already processed, {len(sdf_files)} remaining")
    
    # Create output directory
    Path(output_dir).mkdir(parents=True, exist_ok=True)
    
    # Determine cores
    if num_cores is None:
        num_cores = max(1, cpu_count() - 1)
    
    logger.log(f"Processing parameters:")
    logger.log(f"  - Total files: {len(all_sdf_files)}")
    logger.log(f"  - Already processed: {already_processed}")
    logger.log(f"  - To process: {len(sdf_files)}")
    logger.log(f"  - CPU cores: {num_cores}")
    logger.log(f"  - Output directory: {output_dir}\n")
    
    if len(sdf_files) == 0:
        logger.log("All files already processed!")
        return
    
    success_count = 0
    fail_count = 0
    failed_files = []
    
    # Process files
    logger.log("Starting conversion...\n")
    
    if num_cores > 1:
        # Parallel processing with checkpoint after each batch
        batch_size = 100
        for i in range(0, len(sdf_files), batch_size):
            batch = sdf_files[i:i+batch_size]
            args_list = [(sdf_file, output_dir) for sdf_file in batch]
            
            with Pool(num_cores) as pool:
                results_list = pool.map(convert_single_file_meeko, args_list)
            
            # Process results and checkpoint
            for sdf_file, results in zip(batch, results_list):
                for filename, success, error in results:
                    if success:
                        logger.log(f"✓ [{success_count + fail_count + 1}/{len(sdf_files)}] {filename}")
                        logger.log_success(sdf_file.name, filename)
                        success_count += 1
                        checkpoint_mgr.save(sdf_file.name, True)
                    else:
                        logger.log(f"✗ [{success_count + fail_count + 1}/{len(sdf_files)}] {filename} - {error}")
                        logger.log_error(filename, error)
                        fail_count += 1
                        failed_files.append((filename, error))
                        checkpoint_mgr.save(sdf_file.name, False)
            
            # Progress update
            processed = success_count + fail_count
            logger.log(f"\nProgress: {processed}/{len(sdf_files)} ({processed/len(sdf_files)*100:.1f}%)")
            logger.log(f"Success: {success_count} | Failed: {fail_count}\n")
    else:
        # Single-threaded
        for idx, sdf_file in enumerate(sdf_files, 1):
            results = convert_single_file_meeko((sdf_file, output_dir))
            for filename, success, error in results:
                if success:
                    logger.log(f"✓ [{idx}/{len(sdf_files)}] {filename}")
                    logger.log_success(sdf_file.name, filename)
                    success_count += 1
                    checkpoint_mgr.save(sdf_file.name, True)
                else:
                    logger.log(f"✗ [{idx}/{len(sdf_files)}] {filename} - {error}")
                    logger.log_error(filename, error)
                    fail_count += 1
                    failed_files.append((filename, error))
                    checkpoint_mgr.save(sdf_file.name, False)
    
    # Final summary
    total_processed = already_processed + success_count + fail_count
    total_success = already_processed + success_count  # Assume previously processed were successful
    
    logger.log(f"\n{'='*70}")
    logger.log(f"CONVERSION SUMMARY:")
    logger.log(f"{'='*70}")
    logger.log(f"Total files: {len(all_sdf_files)}")
    logger.log(f"Previously processed: {already_processed}")
    logger.log(f"Newly processed: {success_count + fail_count}")
    logger.log(f"  - Successful: {success_count}")
    logger.log(f"  - Failed: {fail_count}")
    logger.log(f"Overall success rate: {(total_success/len(all_sdf_files)*100):.1f}%")
    logger.log(f"{'='*70}")
    
    # Log failed files
    if failed_files:
        logger.log(f"\nFailed files ({len(failed_files)}):")
        for filename, error in failed_files[:20]:  # Show first 20
            logger.log(f"  - {filename}: {error}")
        if len(failed_files) > 20:
            logger.log(f"  ... and {len(failed_files) - 20} more (see error log for all details)")


def main():
    parser = argparse.ArgumentParser(
        description='Convert SDF files to PDBQT format with checkpointing and separate logs'
    )
    parser.add_argument('-i', '--input', required=True, help='Input directory with SDF files')
    parser.add_argument('-o', '--output', required=True, help='Output directory for PDBQT files')
    parser.add_argument('-c', '--cores', type=int, default=None, help='CPU cores (default: all-1)')
    parser.add_argument('--resume', action='store_true', help='Resume from checkpoint')
    parser.add_argument('--reset', action='store_true', help='Reset checkpoint and start fresh')
    
    args = parser.parse_args()
    
    if not os.path.isdir(args.input):
        print(f"Error: Input directory '{args.input}' does not exist")
        sys.exit(1)
    
    # Setup session name based on input/output dirs
    session_name = f"{Path(args.input).name}_to_{Path(args.output).name}"
    timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
    
    # Create separate log files
    main_log = LOGS_DIR / f"{session_name}_{timestamp}_MAIN.log"
    success_log = LOGS_DIR / f"{session_name}_{timestamp}_SUCCESS.log"
    error_log = LOGS_DIR / f"{session_name}_{timestamp}_ERROR.log"
    checkpoint_file = LOGS_DIR / f"{session_name}_checkpoint.json"
    
    # Reset checkpoint if requested
    if args.reset and checkpoint_file.exists():
        checkpoint_file.unlink()
        print(f"Checkpoint reset: {checkpoint_file}")
    
    # Initialize
    logger = Logger(main_log, success_log, error_log)
    checkpoint_mgr = CheckpointManager(checkpoint_file)
    
    logger.log(f"Configuration:")
    logger.log(f"  Input: {args.input}")
    logger.log(f"  Output: {args.output}")
    logger.log(f"  Main log: {main_log}")
    logger.log(f"  Success log: {success_log}")
    logger.log(f"  Error log: {error_log}")
    logger.log(f"  Checkpoint: {checkpoint_file}")
    logger.log(f"  Resume mode: {'Yes' if args.resume else 'No'}\n")
    
    # Convert
    try:
        Path(args.output).mkdir(parents=True, exist_ok=True)
        convert_with_meeko(args.input, args.output, args.cores, logger, checkpoint_mgr)
    except KeyboardInterrupt:
        logger.log("\n\nConversion interrupted by user. Progress saved to checkpoint.")
        logger.log("Resume with: python sdf_to_pdbqt_converter.py -i INPUT -o OUTPUT --resume")
    except Exception as e:
        logger.log(f"\n\nFATAL ERROR: {e}")
        import traceback
        logger.log(traceback.format_exc())
    finally:
        logger.close()
    
    print(f"\n{'='*70}")
    print(f"Logs saved in: {LOGS_DIR}")
    print(f"  Main log: {main_log.name}")
    print(f"  Success log: {success_log.name}")
    print(f"  Error log: {error_log.name}")
    print(f"  Checkpoint: {checkpoint_file.name}")
    print(f"{'='*70}")


if __name__ == "__main__":
    main()
