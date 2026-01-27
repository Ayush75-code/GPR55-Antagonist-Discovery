# PyMOL Script to Load Best Poses with GPR55 Receptor
# Run in PyMOL: @load_best_poses.pml

# Load receptor
load "d:\Antigravity AI\Major Project\protein_apo_receptor.pdbqt", GPR55
color gray80, GPR55
show cartoon, GPR55

# Load best ligand poses
load "d:\Antigravity AI\Major Project\results\top_hits_analysis\best_poses_pymol\BEST_P0_compound_11569386_3d_-12.3kcal.pdbqt", best_P0
show sticks, best_P0
color cyan, best_P0

load "d:\Antigravity AI\Major Project\results\top_hits_analysis\best_poses_pymol\BEST_P3_compound_3992081_3d_-8.5kcal.pdbqt", best_P3
show sticks, best_P3
color magenta, best_P3

load "d:\Antigravity AI\Major Project\results\top_hits_analysis\best_poses_pymol\BEST_Interface_CHEMBL432162_-9.6kcal.pdbqt", best_Interface
show sticks, best_Interface
color yellow, best_Interface

# Set up view
zoom
bg_color white
ray 1920, 1080
