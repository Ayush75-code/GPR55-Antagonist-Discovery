# =============================================================================
# COMPLETE HTVS SCRIPT FOR GOOGLE COLAB
# =============================================================================
# Copy this entire file into a Colab notebook cell and run it.
# Shows real-time per-ligand progress like the Phase 1 validation.
# =============================================================================

# === CELL 0: INSTALL VINA (Run once per new Colab) ===
"""
from google.colab import drive
drive.mount('/content/drive')

import os
SOFTWARE_DIR = "/content/drive/MyDrive/Colab_Software"
os.makedirs(SOFTWARE_DIR, exist_ok=True)

if not os.path.exists(f"{SOFTWARE_DIR}/vina"):
    print("Installing AutoDock Vina 1.2.5...")
    !wget -q https://github.com/ccsb-scripps/AutoDock-Vina/releases/download/v1.2.5/vina_1.2.5_linux_x86_64 -O {SOFTWARE_DIR}/vina
    !chmod +x {SOFTWARE_DIR}/vina
    print("✓ Vina installed!")
else:
    print("✓ Vina already installed!")

!{SOFTWARE_DIR}/vina --version
"""

# === CELL 1: COMPLETE HTVS SETUP & EXECUTION ===

from google.colab import drive, runtime
import os
import subprocess
import time
import glob
import json
import sys
from datetime import datetime
from concurrent.futures import ThreadPoolExecutor, as_completed
import multiprocessing

# Mount Drive
drive.mount('/content/drive')

# =============================================================================
# CONFIGURATION - MODIFY THESE AS NEEDED
# =============================================================================
BASE_DIR = "/content/drive/MyDrive/Major Project/MAIN DATA/HTVS"
SOFTWARE_DIR = "/content/drive/MyDrive/Colab_Software"
VINA = f"{SOFTWARE_DIR}/vina"

RECEPTOR = f"{BASE_DIR}/protein_apo_receptor.pdbqt"
LIGANDS_DIR = f"{BASE_DIR}/ligands"
RESULTS_DIR = f"{BASE_DIR}/results"
CHECKPOINT_FILE = f"{RESULTS_DIR}/checkpoint.json"

# Docking Parameters
EXHAUSTIVENESS = 32
NUM_MODES = 1
NUM_WORKERS = multiprocessing.cpu_count() - 2  # Leave 2 CPUs for system

# Target Sites (from your Phase 1 validation)
TARGETS = {
    "P0_Orthosteric": {
        "center": (110.02, 111.68, 81.10),
        "size": (22, 22, 22)
    },
    "P3_Allosteric": {
        "center": (105.58, 102.05, 107.35),
        "size": (20, 20, 20)
    },
    "Interface": {
        "center": (112.678, 92.605, 112.545),
        "size": (25, 25, 25)
    }
}

# Hit threshold (from Phase 1 validation)
HIT_THRESHOLD = -8.5

# =============================================================================
# SETUP
# =============================================================================
os.makedirs(RESULTS_DIR, exist_ok=True)

# Fix Vina permissions
if os.path.exists(VINA):
    os.chmod(VINA, 0o755)
else:
    print("❌ ERROR: Vina not found! Run Cell 0 first.")
    raise SystemExit

# =============================================================================
# CHECKPOINT FUNCTIONS (for resume capability)
# =============================================================================
def load_checkpoint():
    """Load checkpoint from disk"""
    if os.path.exists(CHECKPOINT_FILE):
        try:
            with open(CHECKPOINT_FILE, 'r') as f:
                cp = json.load(f)
                # Convert list back to set for efficiency
                if isinstance(cp.get("completed"), list):
                    cp["completed"] = set(cp["completed"])
                else:
                    cp["completed"] = set(cp.get("completed", []))
                return cp
        except:
            pass
    return {"completed": set(), "results": {}, "hits": {}}

def save_checkpoint(cp):
    """Save checkpoint to disk"""
    save_data = {
        "completed": list(cp["completed"]),  # Convert set to list for JSON
        "results": cp["results"],
        "hits": cp["hits"],
        "last_update": datetime.now().isoformat()
    }
    with open(CHECKPOINT_FILE, 'w') as f:
        json.dump(save_data, f)

def is_completed(cp, target, ligand_name):
    """Check if a ligand-target pair is already done"""
    key = f"{target}_{ligand_name}"
    return key in cp["completed"]

def mark_completed(cp, target, ligand_name, affinity):
    """Mark a ligand-target pair as completed"""
    key = f"{target}_{ligand_name}"
    cp["completed"].add(key)
    if affinity is not None:
        cp["results"][key] = affinity
        if affinity < HIT_THRESHOLD:
            cp["hits"][key] = affinity

# =============================================================================
# DOCKING FUNCTION
# =============================================================================
def dock_ligand(ligand_path, target_name, target_config, output_dir):
    """
    Run AutoDock Vina on a single ligand against a target site.
    Returns: (ligand_name, affinity or None)
    """
    ligand_name = os.path.basename(ligand_path).replace('.pdbqt', '')
    output_file = f"{output_dir}/{ligand_name}_out.pdbqt"
    
    cx, cy, cz = target_config["center"]
    sx, sy, sz = target_config["size"]
    
    cmd = [
        VINA,
        "--receptor", RECEPTOR,
        "--ligand", ligand_path,
        "--center_x", str(cx),
        "--center_y", str(cy),
        "--center_z", str(cz),
        "--size_x", str(sx),
        "--size_y", str(sy),
        "--size_z", str(sz),
        "--exhaustiveness", str(EXHAUSTIVENESS),
        "--num_modes", str(NUM_MODES),
        "--out", output_file
    ]
    
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=300)
        
        # Parse affinity from output
        for line in result.stdout.split('\n'):
            if line.strip().startswith('1'):
                parts = line.split()
                if len(parts) >= 2:
                    try:
                        return ligand_name, float(parts[1])
                    except ValueError:
                        pass
                break
        
        return ligand_name, None
        
    except subprocess.TimeoutExpired:
        return ligand_name, None
    except Exception as e:
        return ligand_name, None

# =============================================================================
# MAIN HTVS FUNCTION
# =============================================================================
def run_htvs(target_name=None, terminate_after=False):
    """
    Run High-Throughput Virtual Screening.
    
    Args:
        target_name: Optional, run only this target (e.g., "P0_Orthosteric")
        terminate_after: If True, terminate Colab runtime after completion
    """
    
    # Load checkpoint
    checkpoint = load_checkpoint()
    
    # Get all ligands
    all_ligands = sorted(glob.glob(f"{LIGANDS_DIR}/*.pdbqt"))
    total_ligands = len(all_ligands)
    
    # Determine which targets to run
    if target_name:
        if target_name not in TARGETS:
            print(f"❌ Unknown target: {target_name}")
            print(f"   Available: {list(TARGETS.keys())}")
            return
        targets_to_run = {target_name: TARGETS[target_name]}
    else:
        targets_to_run = TARGETS
    
    # Print header
    print("=" * 80)
    print("     AutoDock Vina - High-Throughput Virtual Screening")
    print("=" * 80)
    print(f"Receptor:        {os.path.basename(RECEPTOR)}")
    print(f"Ligands:         {total_ligands:,}")
    print(f"Targets:         {list(targets_to_run.keys())}")
    print(f"Exhaustiveness:  {EXHAUSTIVENESS}")
    print(f"Modes:           {NUM_MODES}")
    print(f"Workers:         {NUM_WORKERS}")
    print(f"Hit threshold:   {HIT_THRESHOLD} kcal/mol")
    print(f"Resume:          {len(checkpoint['completed']):,} completed")
    print("=" * 80)
    print()
    
    global_start = time.time()
    
    # Process each target
    for tgt_name, tgt_config in targets_to_run.items():
        
        # Create output directory
        target_output_dir = f"{RESULTS_DIR}/{tgt_name}"
        os.makedirs(target_output_dir, exist_ok=True)
        
        # Filter out already completed ligands
        remaining_ligands = []
        for lig_path in all_ligands:
            lig_name = os.path.basename(lig_path).replace('.pdbqt', '')
            if not is_completed(checkpoint, tgt_name, lig_name):
                remaining_ligands.append(lig_path)
        
        pending_count = len(remaining_ligands)
        completed_count = total_ligands - pending_count
        
        print(f"\n🎯 {tgt_name}")
        print(f"   Total: {total_ligands:,} | Done: {completed_count:,} | Remaining: {pending_count:,}")
        
        if pending_count == 0:
            print("   ✅ Already complete!")
            continue
        
        print()
        
        # Tracking variables
        start_time = time.time()
        done_count = 0
        best_affinity = None
        hit_count = 0
        
        # Run docking in parallel
        with ThreadPoolExecutor(max_workers=NUM_WORKERS) as executor:
            
            # Submit all jobs
            futures = {
                executor.submit(dock_ligand, lig, tgt_name, tgt_config, target_output_dir): lig 
                for lig in remaining_ligands
            }
            
            # Process results as they complete
            for future in as_completed(futures):
                ligand_name, affinity = future.result()
                done_count += 1
                
                # Update tracking
                if affinity is not None:
                    if best_affinity is None or affinity < best_affinity:
                        best_affinity = affinity
                    if affinity < HIT_THRESHOLD:
                        hit_count += 1
                
                # Mark as completed
                mark_completed(checkpoint, tgt_name, ligand_name, affinity)
                
                # Calculate stats
                elapsed = time.time() - start_time
                rate = done_count / elapsed if elapsed > 0 else 0
                eta_minutes = (pending_count - done_count) / rate / 60 if rate > 0 else 0
                percent = done_count / pending_count * 100
                
                # Format affinity
                aff_str = f"{affinity:.2f}" if affinity else "FAIL"
                best_str = f"{best_affinity:.2f}" if best_affinity else "N/A"
                
                # Hit indicator
                hit_marker = " 🌟 HIT!" if affinity and affinity < HIT_THRESHOLD else ""
                new_best = " 📈 NEW BEST!" if affinity == best_affinity and done_count > 1 else ""
                
                # Print progress line
                print(f"[{done_count:5d}/{pending_count}] {percent:5.1f}% | "
                      f"✓ {ligand_name[:22]:22s} | "
                      f"Affinity: {aff_str:>7s} | "
                      f"Best: {best_str:>7s} | "
                      f"Hits: {hit_count:4d} | "
                      f"{rate:.1f}/s | "
                      f"ETA: {eta_minutes:.0f}m"
                      f"{hit_marker}{new_best}")
                sys.stdout.flush()
                
                # Save checkpoint periodically
                if done_count % 100 == 0:
                    save_checkpoint(checkpoint)
        
        # Final save for this target
        save_checkpoint(checkpoint)
        
        # Target summary
        elapsed_minutes = (time.time() - start_time) / 60
        print()
        print(f"   ═══════════════════════════════════════════════════════")
        print(f"   ✅ {tgt_name} COMPLETE")
        print(f"   Time: {elapsed_minutes:.1f} min | Best: {best_affinity:.2f} kcal/mol | Hits: {hit_count}")
        print(f"   ═══════════════════════════════════════════════════════")
    
    # Final summary
    total_elapsed = (time.time() - global_start) / 60
    total_completed = len(checkpoint['completed'])
    total_hits = len(checkpoint['hits'])
    
    print()
    print("=" * 80)
    print("                    HTVS COMPLETE!")
    print("=" * 80)
    print(f"Total time:      {total_elapsed:.1f} minutes ({total_elapsed/60:.1f} hours)")
    print(f"Total dockings:  {total_completed:,}")
    print(f"Total hits:      {total_hits:,} (< {HIT_THRESHOLD} kcal/mol)")
    print(f"Results saved:   {RESULTS_DIR}")
    print(f"Checkpoint:      {CHECKPOINT_FILE}")
    print("=" * 80)
    
    # Print top hits
    if checkpoint['hits']:
        print()
        print("🏆 TOP 10 HITS:")
        sorted_hits = sorted(checkpoint['hits'].items(), key=lambda x: x[1])[:10]
        for i, (key, aff) in enumerate(sorted_hits, 1):
            target, ligand = key.split('_', 1)
            print(f"   {i:2d}. {ligand:30s} | {target:15s} | {aff:.2f} kcal/mol")
    
    if terminate_after:
        print()
        print("⏹️ Terminating runtime in 10 seconds...")
        time.sleep(10)
        runtime.unassign()

# =============================================================================
# VERIFICATION FUNCTION
# =============================================================================
def verify_setup():
    """Verify all files and settings before running HTVS"""
    
    print("=" * 80)
    print("HTVS SETUP VERIFICATION")
    print("=" * 80)
    
    all_ok = True
    
    # Check Vina
    print("\n🔧 Vina:")
    if os.path.exists(VINA):
        result = subprocess.run([VINA, "--version"], capture_output=True, text=True)
        version = result.stdout.strip().split('\n')[0] if result.stdout else "Unknown"
        print(f"   ✅ {version}")
    else:
        print(f"   ❌ NOT FOUND at {VINA}")
        all_ok = False
    
    # Check Receptor
    print("\n📄 Receptor:")
    if os.path.exists(RECEPTOR):
        size_kb = os.path.getsize(RECEPTOR) / 1024
        print(f"   ✅ {os.path.basename(RECEPTOR)} ({size_kb:.1f} KB)")
    else:
        print(f"   ❌ NOT FOUND at {RECEPTOR}")
        all_ok = False
    
    # Check Ligands
    print("\n📦 Ligands:")
    if os.path.exists(LIGANDS_DIR):
        ligand_count = len(glob.glob(f"{LIGANDS_DIR}/*.pdbqt"))
        print(f"   ✅ {ligand_count:,} PDBQT files in {LIGANDS_DIR}")
    else:
        print(f"   ❌ NOT FOUND at {LIGANDS_DIR}")
        all_ok = False
    
    # Check Targets
    print("\n🎯 Targets:")
    for name, cfg in TARGETS.items():
        print(f"   ✅ {name}: center={cfg['center']}, size={cfg['size']}")
    
    # System info
    print("\n🖥️ System:")
    print(f"   ✅ CPUs: {multiprocessing.cpu_count()} (using {NUM_WORKERS})")
    
    # Checkpoint
    print("\n💾 Checkpoint:")
    cp = load_checkpoint()
    print(f"   ✅ {len(cp['completed']):,} completed, {len(cp['hits']):,} hits")
    
    # Test single docking
    if all_ok:
        print("\n🧪 Test Docking:")
        test_ligand = glob.glob(f"{LIGANDS_DIR}/*.pdbqt")[0]
        test_name = os.path.basename(test_ligand)
        
        start = time.time()
        name, aff = dock_ligand(test_ligand, "P0_Orthosteric", TARGETS["P0_Orthosteric"], "/tmp")
        elapsed = time.time() - start
        
        if aff:
            print(f"   ✅ {test_name}: {aff:.2f} kcal/mol ({elapsed:.1f}s)")
            
            # Estimate time
            per_ligand = elapsed
            total_jobs = ligand_count * len(TARGETS)
            est_hours = (per_ligand * total_jobs / NUM_WORKERS) / 3600
            print(f"\n⏱️ Time Estimate:")
            print(f"   Per ligand: ~{per_ligand:.1f}s")
            print(f"   Total jobs: {total_jobs:,}")
            print(f"   Estimated:  ~{est_hours:.1f} hours")
        else:
            print(f"   ❌ Test failed!")
            all_ok = False
    
    print()
    print("=" * 80)
    if all_ok:
        print("✅ ALL CHECKS PASSED - Ready to run!")
        print()
        print("To start HTVS, run one of:")
        print("   run_htvs()                          # Run all targets")
        print("   run_htvs('P0_Orthosteric')          # Run single target")
        print("   run_htvs(terminate_after=True)      # Terminate after completion")
    else:
        print("❌ SETUP ISSUES - Fix errors above before running")
    print("=" * 80)


# =============================================================================
# RUN VERIFICATION AUTOMATICALLY
# =============================================================================
if __name__ == "__main__":
    verify_setup()
