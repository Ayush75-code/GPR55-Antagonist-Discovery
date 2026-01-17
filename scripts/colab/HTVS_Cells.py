# =============================================================================
# AUTODOCK-GPU HTVS - L4 Optimized with Resume Support
# =============================================================================
# 
# Features:
#   ✅ 3 batches per session
#   ✅ Resume from interruption (checkpoint system)
#   ✅ Real-time results saved to Google Drive
#   ✅ Progress tracking with ETA
#
# Sessions:
#   Session 1: Batches 01-03
#   Session 2: Batches 04-06
#   Session 3: Batches 07-09
#   Session 4: Batches 10-12
#   Session 5: Batches 13-15
#
# =============================================================================


# %% [markdown]
# # AutoDock-GPU HTVS - GPR55 Allosteric Screening (L4 Optimized)
# 
# **Features:** Resume support, 3 batches/session, real-time saves
# 
# **Targets:** P0 (Orthosteric), P3 (Allosteric), Interface  
# **Exhaustiveness:** 32 | **Runs:** 3/site | **GPU:** L4


# %% [markdown]
# ## CELL 0: One-Time Setup (Install AutoDock-GPU)

# %%
# === CELL 0: ONE-TIME SETUP ===
from google.colab import drive
drive.mount('/content/drive')

import os

SOFTWARE_DIR = "/content/drive/MyDrive/Colab_Software"
AUTODOCK_DIR = f"{SOFTWARE_DIR}/AutoDock-GPU"

os.makedirs(AUTODOCK_DIR, exist_ok=True)

print("=" * 80)
print("INSTALLING AUTODOCK-GPU TO GOOGLE DRIVE")
print("=" * 80)

if os.path.exists(f"{AUTODOCK_DIR}/bin/autodock_gpu_128wi"):
    print("✅ AutoDock-GPU already installed!")
else:
    !apt-get update -qq
    !apt-get install -y -qq build-essential ocl-icd-opencl-dev
    
    %cd /content
    !git clone https://github.com/ccsb-scripps/AutoDock-GPU.git autodock-gpu-src
    %cd /content/autodock-gpu-src
    !make DEVICE=CUDA NUMWI=128
    
    !cp -r /content/autodock-gpu-src/bin {AUTODOCK_DIR}/
    print("✅ AutoDock-GPU installed!")

!ls -la {AUTODOCK_DIR}/bin/ | head -5
print("\n✅ Setup complete! Run Cell 1 next.")


# %% [markdown]
# ## CELL 1: Mount + Setup with Resume Support (Run after every restart)

# %%
# === CELL 1: SETUP WITH RESUME SUPPORT ===
from google.colab import drive, runtime
import os, subprocess, time, glob, json
from datetime import datetime

drive.mount('/content/drive')

# === CONFIGURATION ===
BASE_DIR = "/content/drive/MyDrive/Major Project/MAIN DATA"
SOFTWARE_DIR = "/content/drive/MyDrive/Colab_Software/AutoDock-GPU"
RESULTS_DIR = f"{BASE_DIR}/HTVS_results"
CHECKPOINT_FILE = f"{RESULTS_DIR}/checkpoint.json"
RECEPTOR = f"{BASE_DIR}/protein_apo_receptor.pdbqt"
BATCHES_DIR = f"{BASE_DIR}/ligand_batches"

EXHAUSTIVENESS = 32
NRUN = 3

TARGETS = {
    "P0_Orthosteric": {"center": (110.02, 111.68, 81.10), "size": (22, 22, 22)},
    "P3_Allosteric": {"center": (105.58, 102.05, 107.35), "size": (20, 20, 20)},
    "Interface": {"center": (112.678, 92.605, 112.545), "size": (25, 25, 25)}
}

os.makedirs(RESULTS_DIR, exist_ok=True)
os.environ["PATH"] = f"{SOFTWARE_DIR}/bin:" + os.environ["PATH"]

# === CHECKPOINT FUNCTIONS ===
def load_checkpoint():
    """Load checkpoint from Google Drive"""
    if os.path.exists(CHECKPOINT_FILE):
        with open(CHECKPOINT_FILE, 'r') as f:
            return json.load(f)
    return {"completed": {}, "last_update": None}

def save_checkpoint(checkpoint):
    """Save checkpoint to Google Drive (real-time)"""
    checkpoint["last_update"] = datetime.now().isoformat()
    with open(CHECKPOINT_FILE, 'w') as f:
        json.dump(checkpoint, f, indent=2)

def is_ligand_done(checkpoint, batch, target, ligand):
    """Check if a specific ligand-target combo is already done"""
    key = f"{batch}_{target}_{ligand}"
    return key in checkpoint.get("completed", {})

def mark_ligand_done(checkpoint, batch, target, ligand, affinity):
    """Mark a ligand as completed"""
    key = f"{batch}_{target}_{ligand}"
    checkpoint["completed"][key] = {"affinity": affinity, "time": datetime.now().isoformat()}

# === PROGRESS TRACKING ===
class ProgressTracker:
    def __init__(self, total):
        self.total = total
        self.done = 0
        self.start_time = time.time()
        self.affinities = []
    
    def update(self, affinity=None):
        self.done += 1
        if affinity: self.affinities.append(affinity)
        
        elapsed = time.time() - self.start_time
        rate = self.done / elapsed if elapsed > 0 else 0
        eta = (self.total - self.done) / rate if rate > 0 else 0
        pct = self.done * 100 / self.total
        
        best = min(self.affinities) if self.affinities else 0
        
        if self.done % 50 == 0 or self.done == self.total:
            print(f"   [{self.done:,}/{self.total:,}] {pct:.1f}% | "
                  f"Rate: {rate:.1f}/s | ETA: {eta/60:.1f}min | Best: {best:.2f}")

# === MAIN DOCKING FUNCTION ===
def run_batches(batch_start, batch_end, terminate=True):
    """Run multiple batches with resume support"""
    checkpoint = load_checkpoint()
    
    print("=" * 80)
    print(f"HTVS: BATCHES {batch_start:02d} - {batch_end:02d}")
    print("=" * 80)
    print(f"Resume mode: {len(checkpoint.get('completed', {}))} ligands already done")
    
    total_start = time.time()
    
    for batch_num in range(batch_start, batch_end + 1):
        batch_name = f"batch_{batch_num:02d}"
        batch_dir = f"{BATCHES_DIR}/{batch_name}"
        batch_results = f"{RESULTS_DIR}/{batch_name}"
        
        if not os.path.exists(batch_dir):
            print(f"⚠️ {batch_name} not found, skipping...")
            continue
        
        print(f"\n{'='*60}")
        print(f"📦 {batch_name.upper()}")
        print('='*60)
        
        ligands = sorted(glob.glob(f"{batch_dir}/*.pdbqt"))
        os.makedirs(batch_results, exist_ok=True)
        
        for target_name, cfg in TARGETS.items():
            cx, cy, cz = cfg["center"]
            target_out = f"{batch_results}/{target_name}"
            os.makedirs(target_out, exist_ok=True)
            
            # Count remaining
            remaining = [l for l in ligands if not is_ligand_done(
                checkpoint, batch_name, target_name, os.path.basename(l).replace('.pdbqt', ''))]
            
            print(f"\n🎯 {target_name}: {len(remaining)}/{len(ligands)} remaining")
            
            if len(remaining) == 0:
                print("   ✅ Already complete!")
                continue
            
            tracker = ProgressTracker(len(remaining))
            
            for lig in remaining:
                name = os.path.basename(lig).replace('.pdbqt', '')
                outfile = f"{target_out}/{name}"
                
                # Run AutoDock-GPU
                cmd = f"autodock_gpu_128wi --ffile {RECEPTOR} --lfile {lig} " \
                      f"--nrun {NRUN} --nev 2500000 --resnam {outfile} 2>/dev/null"
                subprocess.run(cmd, shell=True, capture_output=True)
                
                # Parse result
                affinity = None
                dlg_file = f"{outfile}.dlg"
                if os.path.exists(dlg_file):
                    try:
                        with open(dlg_file) as f:
                            for line in f:
                                if 'RANKING' in line and 'kcal' in line.lower():
                                    affinity = float(line.split()[3])
                                    break
                    except: pass
                
                # Update checkpoint (real-time save every 10 ligands)
                mark_ligand_done(checkpoint, batch_name, target_name, name, affinity)
                if tracker.done % 10 == 0:
                    save_checkpoint(checkpoint)
                
                tracker.update(affinity)
            
            # Save at end of target
            save_checkpoint(checkpoint)
        
        print(f"\n✅ {batch_name} complete!")
    
    # Final save
    save_checkpoint(checkpoint)
    
    elapsed = (time.time() - total_start) / 60
    print("\n" + "=" * 80)
    print(f"✅ BATCHES {batch_start:02d}-{batch_end:02d} COMPLETED!")
    print("=" * 80)
    print(f"⏱️ Total time: {elapsed:.1f} minutes")
    print(f"💾 Checkpoint saved: {len(checkpoint['completed']):,} total ligands done")
    print(f"📁 Results: {RESULTS_DIR}")
    
    if terminate:
        print("\n⏹️ Terminating runtime in 5 seconds...")
        time.sleep(5)
        runtime.unassign()

# Verify setup
print("=" * 80)
print("HTVS SETUP")
print("=" * 80)
print(f"✅ Receptor: {os.path.exists(RECEPTOR)}")
batches = sorted([f for f in os.listdir(BATCHES_DIR) if f.startswith('batch_')])
print(f"✅ Batches: {len(batches)}")
!nvidia-smi --query-gpu=name,memory.total --format=csv,noheader

checkpoint = load_checkpoint()
print(f"✅ Checkpoint: {len(checkpoint.get('completed', {}))} ligands already done")
print("\n✅ Setup complete! Run appropriate session cell.")


# %% [markdown]
# ## SESSION CELLS (3 batches each)

# %%
# === SESSION 1: BATCHES 01-03 ===
run_batches(1, 3)

# %%
# === SESSION 2: BATCHES 04-06 ===
run_batches(4, 6)

# %%
# === SESSION 3: BATCHES 07-09 ===
run_batches(7, 9)

# %%
# === SESSION 4: BATCHES 10-12 ===
run_batches(10, 12)

# %%
# === SESSION 5: BATCHES 13-15 ===
run_batches(13, 15)


# %% [markdown]
# ## Check Progress (Run anytime)

# %%
# === CHECK PROGRESS ===
import json, os

CHECKPOINT_FILE = "/content/drive/MyDrive/Major Project/MAIN DATA/HTVS_results/checkpoint.json"

if os.path.exists(CHECKPOINT_FILE):
    with open(CHECKPOINT_FILE) as f:
        cp = json.load(f)
    
    total = len(cp.get('completed', {}))
    expected = 35002 * 3  # 35002 ligands × 3 targets
    pct = total * 100 / expected
    
    print("=" * 60)
    print("HTVS PROGRESS")
    print("=" * 60)
    print(f"Completed: {total:,} / {expected:,} ({pct:.1f}%)")
    print(f"Last update: {cp.get('last_update', 'N/A')}")
    
    # Count per batch
    batch_counts = {}
    for key in cp.get('completed', {}):
        batch = key.split('_')[0] + '_' + key.split('_')[1]
        batch_counts[batch] = batch_counts.get(batch, 0) + 1
    
    print("\nPer batch:")
    for b in sorted(batch_counts.keys()):
        print(f"  {b}: {batch_counts[b]:,}")
else:
    print("No checkpoint found - HTVS not started yet")


# %% [markdown]
# ## Results Aggregation (Run after all complete)

# %%
# === RESULTS AGGREGATION ===
import os, json, pandas as pd

RESULTS_DIR = "/content/drive/MyDrive/Major Project/MAIN DATA/HTVS_results"
CHECKPOINT_FILE = f"{RESULTS_DIR}/checkpoint.json"

print("=" * 60)
print("AGGREGATING RESULTS FROM CHECKPOINT")
print("=" * 60)

with open(CHECKPOINT_FILE) as f:
    cp = json.load(f)

results = []
for key, data in cp.get('completed', {}).items():
    parts = key.split('_')
    batch = f"{parts[0]}_{parts[1]}"
    target = '_'.join(parts[2:-1]) if len(parts) > 3 else parts[2]
    ligand = parts[-1]
    
    results.append({
        'Ligand': ligand,
        'Target': target,
        'Batch': batch,
        'Affinity': data.get('affinity')
    })

df = pd.DataFrame(results)
df = df.dropna(subset=['Affinity'])

print(f"Total results: {len(df):,}")

# Save all results
df.to_csv(f"{RESULTS_DIR}/HTVS_ALL.csv", index=False)
print(f"✅ Saved: HTVS_ALL.csv")

# Hits
hits = df[df['Affinity'] < -8.5].sort_values('Affinity')
hits.to_csv(f"{RESULTS_DIR}/HTVS_HITS.csv", index=False)
print(f"✅ Hits (< -8.5 kcal/mol): {len(hits):,}")

# Summary
print("\n📊 Best per target:")
for t in df['Target'].unique():
    best = df[df['Target']==t]['Affinity'].min()
    mean = df[df['Target']==t]['Affinity'].mean()
    print(f"  {t}: Best={best:.2f}, Mean={mean:.2f}")

print("\n🏆 Top 10 Hits:")
print(hits.head(10).to_string(index=False))
