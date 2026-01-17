#!/usr/bin/env python3
"""
Advanced SDF to PDBQT Converter
Handles multi-fragment molecules by keeping only the largest fragment
For processing problematic molecules that fail in standard converter
"""

import os
import sys
import argparse
from pathlib import Path
from multiprocessing import Pool, cpu_count
from datetime import datetime

# Create logs directory
SCRIPT_DIR = Path(__file__).parent.absolute()
LOGS_DIR = SCRIPT_DIR / "conversion_logs"
LOGS_DIR.mkdir(exist_ok=True)


class Logger:
    """Logger for multiple log files"""
    def __init__(self, main_log, success_log, error_log):
        self.main_log = main_log
        self.success_log = success_log
        self.error_log = error_log
        
        self.main_handle = open(main_log, 'a')
        self.success_handle = open(success_log, 'a')
        self.error_handle = open(error_log, 'a')
        
        timestamp = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
        header = f"=== Advanced SDF to PDBQT Conversion (Multi-Fragment Handler) ===\nSession started at: {timestamp}\n"
        self.main_handle.write(header)
        self.main_handle.flush()
    
    def log(self, message):
        print(message)
        self.main_handle.write(message + '\n')
        self.main_handle.flush()
    
    def log_success(self, filename, output_file, fragment_info=""):
        msg = f"✓ {filename} -> {output_file}"
        if fragment_info:
            msg += f" | {fragment_info}"
        self.success_handle.write(f"{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} | {msg}\n")
        self.success_handle.flush()
    
    def log_error(self, filename, error):
        msg = f"✗ {filename} | Error: {error}"
        self.error_handle.write(f"{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} | {msg}\n")
        self.error_handle.flush()
    
    def close(self):
        footer = f"\nSession finished at: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n{'='*70}\n"
        self.main_handle.write(footer)
        self.main_handle.close()
        self.success_handle.close()
        self.error_handle.close()


def convert_single_file_advanced(args):
    """Convert a single SDF file with advanced fragment handling"""
    sdf_file, output_dir = args
    
    try:
        from rdkit import Chem
        from rdkit.Chem import AllChem, Descriptors
        from meeko import MoleculePreparation, PDBQTWriterLegacy
        
        supplier = Chem.SDMolSupplier(str(sdf_file), removeHs=False, sanitize=True)
        results = []
        
        for idx, mol in enumerate(supplier):
            if mol is None:
                results.append((sdf_file.name, False, "Failed to read molecule", ""))
                continue
            
            try:
                # Get number of fragments
                frags = Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=False)
                num_frags = len(frags)
                
                fragment_info = ""
                
                # If multiple fragments, keep only the largest
                if num_frags > 1:
                    # Find largest fragment by heavy atom count
                    largest_frag = max(frags, key=lambda x: x.GetNumHeavyAtoms())
                    mol = largest_frag
                    fragment_info = f"Multi-fragment ({num_frags} frags) - kept largest"
                
                # Add explicit hydrogens
                mol = Chem.AddHs(mol, addCoords=True)
                
                # Generate 3D coords if missing
                if mol.GetNumConformers() == 0:
                    result = AllChem.EmbedMolecule(mol, randomSeed=42)
                    if result != 0:
                        # Try with different parameters
                        result = AllChem.EmbedMolecule(mol, randomSeed=42, useRandomCoords=True)
                        if result != 0:
                            results.append((sdf_file.name, False, "Failed to generate 3D coordinates", ""))
                            continue
                    AllChem.MMFFOptimizeMolecule(mol, maxIters=200)
                
                # Prepare with Meeko
                preparator = MoleculePreparation()
                mol_prepared_list = preparator.prepare(mol)
                
                if not mol_prepared_list or len(mol_prepared_list) == 0:
                    results.append((sdf_file.name, False, "Meeko returned empty result", ""))
                    continue
                
                mol_prepared = mol_prepared_list[0]
                
                # Output filename
                num_mols = len([m for m in Chem.SDMolSupplier(str(sdf_file), removeHs=False) if m is not None])
                if num_mols > 1:
                    output_file = Path(output_dir) / f"{sdf_file.stem}_{idx+1}.pdbqt"
                else:
                    output_file = Path(output_dir) / f"{sdf_file.stem}.pdbqt"
                
                # Write PDBQT
                pdbqt_output = PDBQTWriterLegacy.write_string(mol_prepared)
                
                if isinstance(pdbqt_output, tuple):
                    pdbqt_string = pdbqt_output[0]
                else:
                    pdbqt_string = pdbqt_output
                
                if not pdbqt_string or len(pdbqt_string) == 0:
                    results.append((sdf_file.name, False, "Empty PDBQT output", ""))
                    continue
                
                with open(output_file, 'w') as f:
                    f.write(pdbqt_string)
                
                results.append((output_file.name, True, None, fragment_info))
                
            except Exception as e:
                results.append((sdf_file.name, False, f"Processing error: {str(e)}", ""))
                continue
        
        if not results:
            results.append((sdf_file.name, False, "No valid molecules in file", ""))
        
        return results
        
    except Exception as e:
        return [(sdf_file.name, False, f"File error: {str(e)}", "")]


def convert_advanced(input_dir, output_dir, num_cores, logger):
    """Convert SDF files with advanced preprocessing"""
    try:
        from rdkit import Chem
        from meeko import MoleculePreparation
    except ImportError as e:
        msg = f"Error: Required packages not installed.\nImport error: {e}"
        logger.log(msg)
        sys.exit(1)
    
    logger.log("Using ADVANCED mode: Multi-fragment handler")
    logger.log("This mode extracts the largest fragment from multi-fragment molecules\n")
    
    # Get all SDF files
    sdf_files = sorted(list(Path(input_dir).glob("*.sdf")))
    
    if not sdf_files:
        logger.log(f"No SDF files found in {input_dir}")
        return
    
    Path(output_dir).mkdir(parents=True, exist_ok=True)
    
    if num_cores is None:
        num_cores = max(1, cpu_count() - 1)
    
    logger.log(f"Processing parameters:")
    logger.log(f"  - Total files: {len(sdf_files)}")
    logger.log(f"  - CPU cores: {num_cores}")
    logger.log(f"  - Output directory: {output_dir}\n")
    
    success_count = 0
    fail_count = 0
    multi_frag_count = 0
    failed_files = []
    
    logger.log("Starting advanced conversion...\n")
    
    if num_cores > 1:
        batch_size = 50
        for i in range(0, len(sdf_files), batch_size):
            batch = sdf_files[i:i+batch_size]
            args_list = [(sdf_file, output_dir) for sdf_file in batch]
            
            with Pool(num_cores) as pool:
                results_list = pool.map(convert_single_file_advanced, args_list)
            
            for sdf_file, results in zip(batch, results_list):
                for filename, success, error, frag_info in results:
                    if success:
                        logger.log(f"✓ [{success_count + fail_count + 1}/{len(sdf_files)}] {filename}")
                        logger.log_success(sdf_file.name, filename, frag_info)
                        success_count += 1
                        if "Multi-fragment" in frag_info:
                            multi_frag_count += 1
                    else:
                        logger.log(f"✗ [{success_count + fail_count + 1}/{len(sdf_files)}] {filename} - {error}")
                        logger.log_error(filename, error)
                        fail_count += 1
                        failed_files.append((filename, error))
            
            processed = success_count + fail_count
            logger.log(f"\nProgress: {processed}/{len(sdf_files)} ({processed/len(sdf_files)*100:.1f}%)")
            logger.log(f"Success: {success_count} | Failed: {fail_count} | Multi-frag handled: {multi_frag_count}\n")
    else:
        for idx, sdf_file in enumerate(sdf_files, 1):
            results = convert_single_file_advanced((sdf_file, output_dir))
            for filename, success, error, frag_info in results:
                if success:
                    logger.log(f"✓ [{idx}/{len(sdf_files)}] {filename}")
                    logger.log_success(sdf_file.name, filename, frag_info)
                    success_count += 1
                    if "Multi-fragment" in frag_info:
                        multi_frag_count += 1
                else:
                    logger.log(f"✗ [{idx}/{len(sdf_files)}] {filename} - {error}")
                    logger.log_error(filename, error)
                    fail_count += 1
                    failed_files.append((filename, error))
    
    # Summary
    logger.log(f"\n{'='*70}")
    logger.log(f"ADVANCED CONVERSION SUMMARY:")
    logger.log(f"{'='*70}")
    logger.log(f"Total files: {len(sdf_files)}")
    logger.log(f"Successful conversions: {success_count}")
    logger.log(f"  - Multi-fragment molecules handled: {multi_frag_count}")
    logger.log(f"Failed conversions: {fail_count}")
    logger.log(f"Success rate: {(success_count/len(sdf_files)*100):.1f}%")
    logger.log(f"{'='*70}")
    
    if failed_files:
        logger.log(f"\nRemaining failed files ({len(failed_files)}):")
        for filename, error in failed_files[:10]:
            logger.log(f"  - {filename}: {error}")
        if len(failed_files) > 10:
            logger.log(f"  ... and {len(failed_files) - 10} more (see error log)")


def main():
    parser = argparse.ArgumentParser(
        description='Advanced SDF to PDBQT converter with multi-fragment handling'
    )
    parser.add_argument('-i', '--input', required=True, help='Input directory with failed SDF files')
    parser.add_argument('-o', '--output', required=True, help='Output directory for PDBQT files')
    parser.add_argument('-c', '--cores', type=int, default=None, help='CPU cores (default: all-1)')
    
    args = parser.parse_args()
    
    if not os.path.isdir(args.input):
        print(f"Error: Input directory '{args.input}' does not exist")
        sys.exit(1)
    
    session_name = f"{Path(args.input).name}_ADVANCED"
    timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
    
    main_log = LOGS_DIR / f"{session_name}_{timestamp}_MAIN.log"
    success_log = LOGS_DIR / f"{session_name}_{timestamp}_SUCCESS.log"
    error_log = LOGS_DIR / f"{session_name}_{timestamp}_ERROR.log"
    
    logger = Logger(main_log, success_log, error_log)
    
    logger.log(f"Configuration:")
    logger.log(f"  Input: {args.input}")
    logger.log(f"  Output: {args.output}")
    logger.log(f"  Main log: {main_log}")
    logger.log(f"  Success log: {success_log}")
    logger.log(f"  Error log: {error_log}\n")
    
    try:
        Path(args.output).mkdir(parents=True, exist_ok=True)
        convert_advanced(args.input, args.output, args.cores, logger)
    except KeyboardInterrupt:
        logger.log("\n\nConversion interrupted by user.")
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
    print(f"{'='*70}")


if __name__ == "__main__":
    main()
