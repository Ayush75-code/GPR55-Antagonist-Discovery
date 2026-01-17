# =============================================================================
# VINA HTVS - FAST FILTERING DOCKING (35,002 → 350 Ligands)
# =============================================================================
# 
# Purpose: Quick screening with reduced parameters to filter 35,002 ligands
#          down to top 350 candidates for final high-quality docking
#
# Parameters:
#   - num_modes = 1 (fastest, only best pose)
#   - exhaustiveness = 16 (reduced for speed)
#   - CPU parallel processing (40 cores)
#
# Targets: P0 (Orthosteric), P3 (Allosteric), Interface
#
# Output: Top 350 ligands ranked by best affinity across all sites
#
# =============================================================================
# Based on MP Phase1validation.ipynb structure
# =============================================================================


# %% [markdown]
# # AutoDock Vina HTVS - Fast Filtering Run
# 
# **Goal:** Screen 35,002 ligands → Select top 350 for final docking
# 
# **Parameters:** num_modes=1, exhaustiveness=16 (optimized for speed)
# 
# **Targets:** P0 (Orthosteric), P3 (Allosteric), Interface


# %%
# === CELL 1: MOUNT GOOGLE DRIVE ===
from google.colab import drive
drive.mount('/content/drive')


# %%
# === CELL 2: VERIFY FILES ===
import os

print("Mounting Google Drive...")
GDRIVE_BASE = "/content/drive/MyDrive/Major Project/MAIN DATA/HTVS"

print(f"\n✓ Working directory: {GDRIVE_BASE}")

# Verify structure
if os.path.exists(GDRIVE_BASE):
    print(f"\nContents of {GDRIVE_BASE}:")
    for item in sorted(os.listdir(GDRIVE_BASE)):
        print(f"  {item}")
else:
    print(f"❌ Directory not found: {GDRIVE_BASE}")
    print("Please create this directory and add:")
    print("  - receptor.pdbqt")
    print("  - ligands/ folder with all 35,002 pdbqt files")


# %%
# === CELL 3: CHECK INPUT FILES ===
import os
import glob

print("Checking for required input files in Google Drive...")
print(f"Location: {GDRIVE_BASE}/\n")

# Check receptor
RECEPTOR = f"{GDRIVE_BASE}/receptor.pdbqt"
if os.path.exists(RECEPTOR):
    size = os.path.getsize(RECEPTOR)
    print(f"✓ receptor.pdbqt               ({size:,} bytes)")
else:
    print("✗ receptor.pdbqt               MISSING!")

# Check ligands folder
LIGANDS_DIR = f"{GDRIVE_BASE}/ligands"
if os.path.exists(LIGANDS_DIR):
    ligand_files = glob.glob(f"{LIGANDS_DIR}/*.pdbqt")
    print(f"✓ ligands/                     ({len(ligand_files):,} files)")
else:
    print("✗ ligands/                     MISSING!")

# System info
print("\n" + "="*80)
print("SYSTEM INFORMATION")
print("="*80)
!lscpu | grep -E "^CPU\(s\)|Model name|Thread|Core|Socket"
print("")
!free -h


# %%
# === CELL 4: INSTALL AUTODOCK VINA ===
print("="*80)
print("INSTALLING AUTODOCK VINA")
print("="*80)

# Install dependencies
!apt-get update -qq
!apt-get install -y -qq wget

# Download and install Vina 1.2.5
!wget -q https://github.com/ccsb-scripps/AutoDock-Vina/releases/download/v1.2.5/vina_1.2.5_linux_x86_64 -O vina
!chmod +x vina
!mv vina /usr/local/bin/

# Install GNU Parallel for better job management
!apt-get install -y -qq parallel bc

# Verify installation
!vina --version

print("✓ AutoDock Vina installed successfully!")


# %%
# === CELL 5: CREATE HTVS SCRIPT WITH FAST FILTERING PARAMETERS ===
import os

# Configuration
GDRIVE_BASE = "/content/drive/MyDrive/Major Project/MAIN DATA/HTVS"
RESULTS_DIR = f"{GDRIVE_BASE}/filtering_results"

os.makedirs(RESULTS_DIR, exist_ok=True)

htvs_script = f"""#!/bin/bash

# HTVS Fast Filtering Script - AutoDock Vina
# Parameters: num_modes=1, exhaustiveness=16
# Results saved to: {RESULTS_DIR}/

# ============ CONFIGURATION ============
WORK_DIR="{GDRIVE_BASE}"
RESULTS_DIR="{RESULTS_DIR}"
cd "$WORK_DIR" || exit 1

RECEPTOR="receptor.pdbqt"
LIGANDS_DIR="ligands"
MAX_PARALLEL_JOBS=40

# Fast filtering parameters
NUM_MODES=1
EXHAUSTIVENESS=16

# Target sites (3 main binding sites)
declare -a TARGETS=(
    "P0:center_x=110.02:center_y=111.68:center_z=81.10:size_x=22:size_y=22:size_z=22:Orthosteric"
    "P3:center_x=105.58:center_y=102.05:center_z=107.35:size_x=20:size_y=20:size_z=20:Allosteric"
    "Interface:center_x=112.678:center_y=92.605:center_z=112.545:size_x=25:size_y=25:size_z=25:Interface"
)

# ============ SETUP ============
mkdir -p "$RESULTS_DIR/{{raw_outputs,logs,summary}}"
JOBLIST_FILE="$RESULTS_DIR/joblist.txt"
PROGRESS_FILE="$RESULTS_DIR/progress.count"
CHECKPOINT_FILE="$RESULTS_DIR/checkpoint.txt"

TIMESTAMP=$(date '+%Y%m%d_%H%M%S')
MASTER_LOG="$RESULTS_DIR/logs/master_log_${{TIMESTAMP}}.log"
exec > >(tee -a "$MASTER_LOG") 2>&1

echo "Working directory: $WORK_DIR"
echo "Results directory: $RESULTS_DIR"
echo "Master log file: $MASTER_LOG"
echo "Started at: $(date)"
echo ""

# ============ COUNT LIGANDS ============
TOTAL_LIGANDS=$(ls -1 "$LIGANDS_DIR"/*.pdbqt 2>/dev/null | wc -l)
TOTAL_JOBS=$((TOTAL_LIGANDS * ${{#TARGETS[@]}}))

echo "========================================================================"
echo "     HTVS FAST FILTERING - $TOTAL_LIGANDS Ligands to Top 350"
echo "========================================================================"
echo "Receptor:           $RECEPTOR"
echo "Ligands:            $TOTAL_LIGANDS"
echo "Targets:            ${{#TARGETS[@]}} (P0, P3, Interface)"
echo "Total docking jobs: $TOTAL_JOBS"
echo "Parameters:         num_modes=$NUM_MODES, exhaustiveness=$EXHAUSTIVENESS"
echo "Parallel jobs:      $MAX_PARALLEL_JOBS / $(nproc) cores"
echo "========================================================================"
echo ""

# Validate receptor
if [ ! -f "$RECEPTOR" ]; then
    echo "ERROR: Receptor file not found: $RECEPTOR"
    exit 1
fi
echo "✓ Receptor file validated"

# ============ FUNCTIONS ============

check_existing_results() {{
    echo "Checking for existing results..."
    local completed=0
    
    if [ -f "$CHECKPOINT_FILE" ]; then
        completed=$(wc -l < "$CHECKPOINT_FILE" 2>/dev/null || echo 0)
    fi
    
    echo "Found $completed completed jobs"
    EXISTING_COMPLETED=$completed
}}

generate_joblist() {{
    echo "Generating job list (skipping completed runs)..."
    > "$JOBLIST_FILE"
    
    local skipped=0
    local added=0
    
    for ligand_file in "$LIGANDS_DIR"/*.pdbqt; do
        ligand_name=$(basename "$ligand_file" .pdbqt)
        
        for target in "${{TARGETS[@]}}"; do
            IFS=':' read -r target_id cx cy cz sx sy sz target_name <<< "$target"
            
            # Extract values
            center_x=${{cx#center_x=}}
            center_y=${{cy#center_y=}}
            center_z=${{cz#center_z=}}
            size_x=${{sx#size_x=}}
            size_y=${{sy#size_y=}}
            size_z=${{sz#size_z=}}
            
            local output_file="$RESULTS_DIR/raw_outputs/${{ligand_name}}_${{target_id}}.pdbqt"
            
            if [ -f "$output_file" ] && [ -s "$output_file" ]; then
                ((skipped++))
            else
                echo "$ligand_file|$ligand_name|$target_id|$center_x|$center_y|$center_z|$size_x|$size_y|$size_z|$target_name" >> "$JOBLIST_FILE"
                ((added++))
            fi
        done
    done
    
    echo "Skipped $skipped completed jobs"
    echo "Added $added jobs to process"
}}

run_vina_job() {{
    local job_line=$1
    IFS='|' read -r ligand_file ligand_name target_id center_x center_y center_z size_x size_y size_z target_name <<< "$job_line"
    
    local output_file="$RESULTS_DIR/raw_outputs/${{ligand_name}}_${{target_id}}.pdbqt"
    local log_file="$RESULTS_DIR/raw_outputs/${{ligand_name}}_${{target_id}}.log"
    
    # Run Vina with fast parameters
    if vina --receptor "$RECEPTOR" \\
         --ligand "$ligand_file" \\
         --center_x $center_x \\
         --center_y $center_y \\
         --center_z $center_z \\
         --size_x $size_x \\
         --size_y $size_y \\
         --size_z $size_z \\
         --exhaustiveness $EXHAUSTIVENESS \\
         --num_modes $NUM_MODES \\
         --out "$output_file" \\
         > "$log_file" 2>&1; then
        
        # Parse best affinity
        local best_affinity=$(grep "^   1 " "$log_file" 2>/dev/null | awk '{{print $2}}')
        
        # Record result
        echo "$ligand_name|$target_id|$best_affinity" >> "$RESULTS_DIR/all_results.csv"
        echo "$ligand_name|$target_id" >> "$CHECKPOINT_FILE"
        
        # Update progress
        echo "1" >> "$PROGRESS_FILE"
        local completed=$(wc -l < "$PROGRESS_FILE" 2>/dev/null || echo 0)
        local percent=$(awk "BEGIN {{printf \\"%.1f\\", ($completed/$JOBS_TO_RUN)*100}}")
        
        # Print progress every 100 jobs
        if (( completed % 100 == 0 )); then
            printf "[%6d/%d] %5.1f%% | ✓ %s | %s | Affinity: %s kcal/mol\\n" \\
                   "$completed" "$JOBS_TO_RUN" "$percent" \\
                   "$ligand_name" "$target_id" "$best_affinity"
        fi
    else
        echo "1" >> "$PROGRESS_FILE"
        printf "✗ FAILED: %s %s\\n" "$ligand_name" "$target_id"
    fi
}}

export -f run_vina_job
export RECEPTOR RESULTS_DIR EXHAUSTIVENESS NUM_MODES PROGRESS_FILE JOBS_TO_RUN CHECKPOINT_FILE

# ============ MAIN EXECUTION ============

check_existing_results
generate_joblist

JOBS_TO_RUN=$(wc -l < "$JOBLIST_FILE")
export JOBS_TO_RUN

if [ "$JOBS_TO_RUN" -eq 0 ]; then
    echo "========================================================================"
    echo "✓ All jobs are already completed!"
    echo "========================================================================"
else
    > "$PROGRESS_FILE"
    
    # Initialize results file
    if [ ! -f "$RESULTS_DIR/all_results.csv" ]; then
        echo "Ligand|Target|Affinity" > "$RESULTS_DIR/all_results.csv"
    fi
    
    echo ""
    echo "Starting parallel docking runs..."
    echo "Jobs to process: $JOBS_TO_RUN (resuming from $EXISTING_COMPLETED completed)"
    echo "Using GNU Parallel with $MAX_PARALLEL_JOBS parallel jobs"
    echo "========================================================================"
    echo ""
    
    START_TIME=$(date +%s)
    
    # Run parallel docking
    cat "$JOBLIST_FILE" | parallel --line-buffer -j $MAX_PARALLEL_JOBS run_vina_job {{}}
    
    END_TIME=$(date +%s)
    ELAPSED=$((END_TIME - START_TIME))
    ELAPSED_MIN=$((ELAPSED / 60))
    
    echo ""
    echo "========================================================================"
    echo "✓ HTVS FILTERING COMPLETED!"
    echo "========================================================================"
    echo "Time elapsed: $ELAPSED_MIN minutes ($ELAPSED seconds)"
    echo "Results saved to: $RESULTS_DIR/"
    echo "========================================================================"
fi

echo ""
echo "✓ Docking complete! Run the analysis cell to select top 350."
"""

with open('/tmp/run_htvs_filtering.sh', 'w') as f:
    f.write(htvs_script)

!chmod +x /tmp/run_htvs_filtering.sh

print("✓ HTVS filtering script created!")
print(f"✓ Results will be saved to: {RESULTS_DIR}/")


# %%
# === CELL 6: RUN HTVS FILTERING ===
import time

print("="*80)
print("STARTING HTVS FAST FILTERING")
print("="*80)
print("")
print("Configuration:")
print("  - 35,002 ligands × 3 targets = ~105,006 total dockings")
print("  - num_modes = 1 (fastest)")
print("  - exhaustiveness = 16 (fast filtering)")
print("  - 40 parallel jobs on available cores")
print(f"  - Results saving to: {RESULTS_DIR}/")
print("")
print("⚠️  IMPORTANT: Keep this tab active to prevent disconnection!")
print("   If disconnected, just re-run this cell - auto-resume will continue.")
print("")
print("="*80)
print("")

start_time = time.time()

!/tmp/run_htvs_filtering.sh

end_time = time.time()
elapsed_minutes = (end_time - start_time) / 60

print("")
print("="*80)
print("HTVS FILTERING COMPLETED!")
print("="*80)
print(f"Total time: {elapsed_minutes:.1f} minutes")
print(f"Results saved to: {RESULTS_DIR}/")
print("="*80)


# %%
# === CELL 7: ANALYZE RESULTS & SELECT TOP 350 ===
import pandas as pd
import os

RESULTS_DIR = "/content/drive/MyDrive/Major Project/MAIN DATA/HTVS/filtering_results"

print("="*80)
print("ANALYZING RESULTS - SELECTING TOP 350 LIGANDS")
print("="*80)

# Load all results
results_file = f"{RESULTS_DIR}/all_results.csv"
df = pd.read_csv(results_file, sep='|')

# Remove header row if duplicated
df = df[df['Ligand'] != 'Ligand']
df['Affinity'] = pd.to_numeric(df['Affinity'], errors='coerce')
df = df.dropna(subset=['Affinity'])

print(f"\n📊 Total valid results: {len(df):,}")
print(f"📊 Unique ligands: {df['Ligand'].nunique():,}")

# Summary per target
print("\n📈 Summary per Target:")
print("-" * 50)
for target in ['P0', 'P3', 'Interface']:
    target_df = df[df['Target'] == target]
    if len(target_df) > 0:
        print(f"{target}:")
        print(f"   Count: {len(target_df):,}")
        print(f"   Best: {target_df['Affinity'].min():.2f} kcal/mol")
        print(f"   Mean: {target_df['Affinity'].mean():.2f} kcal/mol")
        print(f"   Worst: {target_df['Affinity'].max():.2f} kcal/mol")

# Get best affinity per ligand (across all targets)
df_best = df.loc[df.groupby('Ligand')['Affinity'].idxmin()]
df_best = df_best.sort_values('Affinity')

print("\n" + "="*80)
print("TOP 350 LIGANDS (Best Affinity Across All Sites)")
print("="*80)

# Select top 350
TOP_N = 350
top_350 = df_best.head(TOP_N).copy()
top_350.reset_index(drop=True, inplace=True)
top_350['Rank'] = top_350.index + 1

print(f"\n📋 Selected {len(top_350)} ligands for final docking")
print(f"   Affinity range: {top_350['Affinity'].min():.2f} to {top_350['Affinity'].max():.2f} kcal/mol")

# Save top 350 list
top_350_file = f"{RESULTS_DIR}/top_350_ligands.csv"
top_350.to_csv(top_350_file, index=False)
print(f"\n💾 Saved: {top_350_file}")

# Save all results sorted
all_sorted = df_best.copy()
all_sorted.reset_index(drop=True, inplace=True)
all_sorted['Rank'] = all_sorted.index + 1
all_sorted.to_csv(f"{RESULTS_DIR}/all_ligands_ranked.csv", index=False)
print(f"💾 Saved: {RESULTS_DIR}/all_ligands_ranked.csv")

# Show top 20
print("\n🏆 TOP 20 HITS:")
print("-" * 60)
print(top_350[['Rank', 'Ligand', 'Target', 'Affinity']].head(20).to_string(index=False))

# Distribution analysis
print("\n📊 Affinity Distribution of Top 350:")
print(f"   < -10.0 kcal/mol: {len(top_350[top_350['Affinity'] < -10.0]):,}")
print(f"   -10.0 to -9.0: {len(top_350[(top_350['Affinity'] >= -10.0) & (top_350['Affinity'] < -9.0)]):,}")
print(f"   -9.0 to -8.5: {len(top_350[(top_350['Affinity'] >= -9.0) & (top_350['Affinity'] < -8.5)]):,}")
print(f"   -8.5 to -8.0: {len(top_350[(top_350['Affinity'] >= -8.5) & (top_350['Affinity'] < -8.0)]):,}")
print(f"   > -8.0 kcal/mol: {len(top_350[top_350['Affinity'] >= -8.0]):,}")


# %%
# === CELL 8: EXPORT TOP 350 LIGAND FILES ===
import shutil

LIGANDS_DIR = "/content/drive/MyDrive/Major Project/MAIN DATA/HTVS/ligands"
RESULTS_DIR = "/content/drive/MyDrive/Major Project/MAIN DATA/HTVS/filtering_results"

print("="*80)
print("EXPORTING TOP 350 LIGAND FILES")
print("="*80)

# Load top 350 list
top_350_file = f"{RESULTS_DIR}/top_350_ligands.csv"
top_350 = pd.read_csv(top_350_file)

# Create output directory for top ligands
TOP_LIGANDS_DIR = "/content/drive/MyDrive/Major Project/MAIN DATA/HTVS/top_350_ligands"
os.makedirs(TOP_LIGANDS_DIR, exist_ok=True)

print(f"\n📂 Source: {LIGANDS_DIR}")
print(f"📂 Destination: {TOP_LIGANDS_DIR}")

# Copy top 350 ligand files
copied = 0
missing = 0

for idx, row in top_350.iterrows():
    ligand_name = row['Ligand']
    src_file = f"{LIGANDS_DIR}/{ligand_name}.pdbqt"
    dst_file = f"{TOP_LIGANDS_DIR}/{ligand_name}.pdbqt"
    
    if os.path.exists(src_file):
        shutil.copy2(src_file, dst_file)
        copied += 1
    else:
        print(f"⚠️ Not found: {ligand_name}.pdbqt")
        missing += 1
    
    # Progress
    if (idx + 1) % 50 == 0:
        print(f"   Copied {idx + 1}/{len(top_350)}...")

print(f"\n✅ Copied: {copied} ligand files")
if missing > 0:
    print(f"⚠️ Missing: {missing} files")

print(f"\n" + "="*80)
print("FILTERING COMPLETE!")
print("="*80)
print(f"📁 Top 350 ligands exported to: {TOP_LIGANDS_DIR}")
print(f"📊 Ready for final high-quality docking with:")
print(f"   - Higher exhaustiveness (32)")
print(f"   - More modes (9)")
print(f"   - Full pose analysis")


# %%
# === CELL 9: CHECK PROGRESS (Run Anytime) ===
import os
import glob

RESULTS_DIR = "/content/drive/MyDrive/Major Project/MAIN DATA/HTVS/filtering_results"
LIGANDS_DIR = "/content/drive/MyDrive/Major Project/MAIN DATA/HTVS/ligands"

print("="*60)
print("HTVS FILTERING PROGRESS")
print("="*60)

checkpoint_file = f"{RESULTS_DIR}/checkpoint.txt"
if os.path.exists(checkpoint_file):
    total_done = sum(1 for _ in open(checkpoint_file))
    
    total_ligands = len(glob.glob(f"{LIGANDS_DIR}/*.pdbqt"))
    total_expected = total_ligands * 3  # 3 targets
    pct = total_done * 100 / total_expected if total_expected > 0 else 0
    
    print(f"\n📊 Progress: {total_done:,} / {total_expected:,} ({pct:.1f}%)")
    print(f"📋 Ligands: {total_ligands:,}")
    print(f"🎯 Targets: 3 (P0, P3, Interface)")
    
    # Count per target
    if os.path.exists(f"{RESULTS_DIR}/all_results.csv"):
        import pandas as pd
        try:
            df = pd.read_csv(f"{RESULTS_DIR}/all_results.csv", sep='|')
            df = df[df['Ligand'] != 'Ligand']
            for target in ['P0', 'P3', 'Interface']:
                count = len(df[df['Target'] == target])
                pct_t = count * 100 / total_ligands if total_ligands > 0 else 0
                print(f"   {target}: {count:,} / {total_ligands:,} ({pct_t:.1f}%)")
        except:
            pass
else:
    print("❌ No checkpoint found - HTVS not started yet")
