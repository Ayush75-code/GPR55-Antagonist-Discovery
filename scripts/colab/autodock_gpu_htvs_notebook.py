# =============================================================================
# AutoDock-GPU HTVS Notebook
# =============================================================================
# 
# This file contains all cells for the Colab notebook.
# Copy each cell into your Colab notebook.
#
# Workflow:
#   Cell 0: One-time setup (install AutoDock-GPU to Google Drive)
#   Cell 1: Mount drive + Setup (run after every restart)
#   Cell 2-16: Batch 01-15 (one per session, ends with runtime termination)
#
# =============================================================================

# =============================================================================
# CELL 0: ONE-TIME SETUP - Install AutoDock-GPU to Google Drive
# =============================================================================
# Run this ONLY ONCE. It installs AutoDock-GPU to your Google Drive
# so you don't need to reinstall every session.
# =============================================================================

"""
# Mount Google Drive
from google.colab import drive
drive.mount('/content/drive')

import os
import subprocess

# Create software directory in Google Drive
SOFTWARE_DIR = "/content/drive/MyDrive/Colab_Software"
AUTODOCK_DIR = f"{SOFTWARE_DIR}/AutoDock-GPU"

os.makedirs(SOFTWARE_DIR, exist_ok=True)
os.makedirs(AUTODOCK_DIR, exist_ok=True)

print("=" * 80)
print("INSTALLING AUTODOCK-GPU TO GOOGLE DRIVE")
print("=" * 80)

# Check if already installed
if os.path.exists(f"{AUTODOCK_DIR}/bin/autodock_gpu_128wi"):
    print("✅ AutoDock-GPU already installed!")
    print(f"   Location: {AUTODOCK_DIR}")
else:
    print("Installing AutoDock-GPU...")
    
    # Install dependencies
    !apt-get update -qq
    !apt-get install -y -qq build-essential ocl-icd-opencl-dev
    
    # Clone and build AutoDock-GPU
    os.chdir("/content")
    !git clone https://github.com/ccsb-scripps/AutoDock-GPU.git autodock-gpu-src
    os.chdir("/content/autodock-gpu-src")
    
    # Build with CUDA support
    !make DEVICE=CUDA NUMWI=128
    
    # Copy to Google Drive
    !cp -r /content/autodock-gpu-src/bin {AUTODOCK_DIR}/
    !cp -r /content/autodock-gpu-src/input {AUTODOCK_DIR}/
    
    print("✅ AutoDock-GPU installed to Google Drive!")

# Verify installation
print("\\n--- Verifying Installation ---")
!ls -la {AUTODOCK_DIR}/bin/

print("\\n✅ Setup complete! You can now run Cell 1 after any restart.")
"""


# =============================================================================
# CELL 1: MOUNT DRIVE + SETUP (Run after every restart)
# =============================================================================
# Run this cell every time you restart the runtime.
# It sets up paths and creates the docking configuration.
# =============================================================================

"""
from google.colab import drive
import os
import subprocess
import time

# Mount Google Drive
print("Mounting Google Drive...")
drive.mount('/content/drive')

# === CONFIGURATION ===
BASE_DIR = "/content/drive/MyDrive/Major Project/MAIN DATA"
SOFTWARE_DIR = "/content/drive/MyDrive/Colab_Software/AutoDock-GPU"
RESULTS_DIR = f"{BASE_DIR}/HTVS_results"

RECEPTOR = f"{BASE_DIR}/protein_apo_receptor.pdbqt"
BATCHES_DIR = f"{BASE_DIR}/ligand_batches"

# Docking parameters
EXHAUSTIVENESS = 32
NRUN = 3  # 3 runs per site

# Target sites
TARGETS = {
    "P0_Orthosteric": {"center": (110.02, 111.68, 81.10), "size": (22, 22, 22)},
    "P3_Allosteric": {"center": (105.58, 102.05, 107.35), "size": (20, 20, 20)},
    "Interface": {"center": (112.678, 92.605, 112.545), "size": (25, 25, 25)}
}

# Create results directory
os.makedirs(RESULTS_DIR, exist_ok=True)

# Add AutoDock-GPU to PATH
os.environ["PATH"] = f"{SOFTWARE_DIR}/bin:" + os.environ["PATH"]

# Verify setup
print("\\n" + "=" * 80)
print("HTVS SETUP VERIFICATION")
print("=" * 80)

print(f"\\n📂 Base Directory: {BASE_DIR}")
print(f"📂 Results Directory: {RESULTS_DIR}")
print(f"🔬 Receptor: {os.path.basename(RECEPTOR)}")

# Check receptor
if os.path.exists(RECEPTOR):
    print("✅ Receptor file found")
else:
    print("❌ ERROR: Receptor file not found!")

# Check batches
batch_folders = sorted([f for f in os.listdir(BATCHES_DIR) if f.startswith("batch_")])
print(f"\\n📦 Found {len(batch_folders)} batch folders:")
for bf in batch_folders:
    count = len([f for f in os.listdir(f"{BATCHES_DIR}/{bf}") if f.endswith('.pdbqt')])
    print(f"   {bf}: {count:,} ligands")

# Check AutoDock-GPU
print("\\n🔧 AutoDock-GPU:")
!which autodock_gpu_128wi || echo "❌ AutoDock-GPU not found in PATH"

# Check GPU
print("\\n🖥️ GPU Status:")
!nvidia-smi --query-gpu=name,memory.total,memory.free --format=csv

print("\\n" + "=" * 80)
print("✅ Setup complete! Run the batch cell you need.")
print("=" * 80)

# Define the docking function for use in batch cells
def run_batch_htvs(batch_num):
    '''Run HTVS for a single batch on all 3 targets'''
    from google.colab import runtime
    import glob
    import shutil
    
    batch_name = f"batch_{batch_num:02d}"
    batch_dir = f"{BATCHES_DIR}/{batch_name}"
    batch_results = f"{RESULTS_DIR}/{batch_name}"
    
    print("=" * 80)
    print(f"HTVS: {batch_name.upper()}")
    print("=" * 80)
    
    if not os.path.exists(batch_dir):
        print(f"❌ ERROR: Batch directory not found: {batch_dir}")
        return
    
    # Get ligands
    ligands = sorted(glob.glob(f"{batch_dir}/*.pdbqt"))
    total_ligands = len(ligands)
    print(f"📋 Ligands to dock: {total_ligands:,}")
    print(f"🎯 Targets: P0, P3, Interface")
    print(f"⚙️ Exhaustiveness: {EXHAUSTIVENESS}")
    print(f"🔄 Runs per site: {NRUN}")
    print(f"📊 Total docking jobs: {total_ligands * 3:,}")
    print("=" * 80)
    
    os.makedirs(batch_results, exist_ok=True)
    start_time = time.time()
    
    # Process each target
    for target_name, target_config in TARGETS.items():
        cx, cy, cz = target_config["center"]
        sx, sy, sz = target_config["size"]
        
        target_results = f"{batch_results}/{target_name}"
        os.makedirs(target_results, exist_ok=True)
        
        print(f"\\n🎯 Target: {target_name}")
        print(f"   Center: ({cx}, {cy}, {cz})")
        print(f"   Size: ({sx}, {sy}, {sz})")
        
        # Create grid parameter file
        gpf_content = f'''npts {sx} {sy} {sz}
gridfld protein.maps.fld
spacing 0.375
receptor_types A C HD N NA OA SA
ligand_types A C F I NA OA N SA S Cl HD Br
gridcenter {cx} {cy} {cz}
smooth 0.5
map protein.A.map
map protein.C.map
map protein.HD.map
map protein.N.map
map protein.NA.map
map protein.OA.map
map protein.SA.map
elecmap protein.e.map
dsolvmap protein.d.map
dielectric -0.1465
'''
        
        # Process ligands
        for i, ligand_path in enumerate(ligands, 1):
            ligand_name = os.path.basename(ligand_path).replace('.pdbqt', '')
            output_prefix = f"{target_results}/{ligand_name}"
            
            # Run AutoDock-GPU
            cmd = f'''autodock_gpu_128wi \\
                --ffile {RECEPTOR} \\
                --lfile {ligand_path} \\
                --nrun {NRUN} \\
                --nev 2500000 \\
                --resnam {output_prefix}'''
            
            result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
            
            # Progress
            if i % 100 == 0 or i == total_ligands:
                elapsed = time.time() - start_time
                rate = i / elapsed if elapsed > 0 else 0
                eta = (total_ligands * 3 - i) / rate if rate > 0 else 0
                print(f"   [{i:,}/{total_ligands:,}] {i*100/total_ligands:.1f}% | Rate: {rate:.1f}/s | ETA: {eta/60:.1f}min")
    
    # Summary
    elapsed = time.time() - start_time
    print("\\n" + "=" * 80)
    print(f"✅ {batch_name.upper()} COMPLETED!")
    print("=" * 80)
    print(f"⏱️ Total time: {elapsed/60:.1f} minutes")
    print(f"📁 Results saved to: {batch_results}")
    print("=" * 80)
    
    # Terminate runtime
    print("\\n🔄 Terminating runtime to save compute units...")
    print("   Re-run Cell 1 and the next batch cell after reconnecting.")
    time.sleep(5)
    runtime.unassign()
"""


# =============================================================================
# CELL 2: BATCH 01
# =============================================================================
"""
run_batch_htvs(1)
"""


# =============================================================================
# CELL 3: BATCH 02
# =============================================================================
"""
run_batch_htvs(2)
"""


# =============================================================================
# CELL 4: BATCH 03
# =============================================================================
"""
run_batch_htvs(3)
"""


# =============================================================================
# CELL 5: BATCH 04
# =============================================================================
"""
run_batch_htvs(4)
"""


# =============================================================================
# CELL 6: BATCH 05
# =============================================================================
"""
run_batch_htvs(5)
"""


# =============================================================================
# CELL 7: BATCH 06
# =============================================================================
"""
run_batch_htvs(6)
"""


# =============================================================================
# CELL 8: BATCH 07
# =============================================================================
"""
run_batch_htvs(7)
"""


# =============================================================================
# CELL 9: BATCH 08
# =============================================================================
"""
run_batch_htvs(8)
"""


# =============================================================================
# CELL 10: BATCH 09
# =============================================================================
"""
run_batch_htvs(9)
"""


# =============================================================================
# CELL 11: BATCH 10
# =============================================================================
"""
run_batch_htvs(10)
"""


# =============================================================================
# CELL 12: BATCH 11
# =============================================================================
"""
run_batch_htvs(11)
"""


# =============================================================================
# CELL 13: BATCH 12
# =============================================================================
"""
run_batch_htvs(12)
"""


# =============================================================================
# CELL 14: BATCH 13
# =============================================================================
"""
run_batch_htvs(13)
"""


# =============================================================================
# CELL 15: BATCH 14
# =============================================================================
"""
run_batch_htvs(14)
"""


# =============================================================================
# CELL 16: BATCH 15 (FINAL)
# =============================================================================
"""
run_batch_htvs(15)
"""


# =============================================================================
# CELL 17: RESULTS AGGREGATION (Run after all batches complete)
# =============================================================================
"""
import os
import glob
import pandas as pd

print("=" * 80)
print("AGGREGATING HTVS RESULTS")
print("=" * 80)

RESULTS_DIR = "/content/drive/MyDrive/Major Project/MAIN DATA/HTVS_results"

all_results = []

# Parse each batch
for batch_num in range(1, 16):
    batch_name = f"batch_{batch_num:02d}"
    batch_dir = f"{RESULTS_DIR}/{batch_name}"
    
    if not os.path.exists(batch_dir):
        print(f"⚠️ {batch_name} not found, skipping...")
        continue
    
    print(f"Processing {batch_name}...")
    
    for target in ["P0_Orthosteric", "P3_Allosteric", "Interface"]:
        target_dir = f"{batch_dir}/{target}"
        if not os.path.exists(target_dir):
            continue
        
        dlg_files = glob.glob(f"{target_dir}/*.dlg")
        
        for dlg_file in dlg_files:
            ligand_name = os.path.basename(dlg_file).replace('.dlg', '')
            
            try:
                with open(dlg_file, 'r') as f:
                    content = f.read()
                
                # Parse best binding energy
                for line in content.split('\\n'):
                    if 'RANKING' in line and 'binding' in line.lower():
                        parts = line.split()
                        energy = float(parts[3])
                        all_results.append({
                            'Ligand': ligand_name,
                            'Target': target,
                            'Batch': batch_name,
                            'Affinity_kcal_mol': energy
                        })
                        break
            except:
                continue

# Create DataFrame
df = pd.DataFrame(all_results)

if len(df) > 0:
    # Best affinity per ligand-target
    df_best = df.loc[df.groupby(['Ligand', 'Target'])['Affinity_kcal_mol'].idxmin()]
    
    # Save results
    df_best.to_csv(f"{RESULTS_DIR}/HTVS_ALL_RESULTS.csv", index=False)
    
    # Summary stats
    print("\\n" + "=" * 80)
    print("HTVS SUMMARY")
    print("=" * 80)
    print(f"\\nTotal ligands processed: {df_best['Ligand'].nunique():,}")
    print(f"Total results: {len(df_best):,}")
    
    print("\\n📊 Best affinities per target:")
    for target in ["P0_Orthosteric", "P3_Allosteric", "Interface"]:
        target_data = df_best[df_best['Target'] == target]
        if len(target_data) > 0:
            best = target_data['Affinity_kcal_mol'].min()
            mean = target_data['Affinity_kcal_mol'].mean()
            print(f"   {target}: Best={best:.2f}, Mean={mean:.2f} kcal/mol")
    
    # Top hits (below threshold)
    threshold = -8.5
    hits = df_best[df_best['Affinity_kcal_mol'] < threshold]
    print(f"\\n🎯 Hits (< {threshold} kcal/mol): {len(hits):,}")
    
    if len(hits) > 0:
        hits.to_csv(f"{RESULTS_DIR}/HTVS_HITS.csv", index=False)
        print(f"   Saved to: {RESULTS_DIR}/HTVS_HITS.csv")
    
    print("\\n✅ Results aggregation complete!")
else:
    print("❌ No results found!")
"""
