# ============================================================
# Figure: HTVS Target Sites (P0, P3, Interface)
# Run in PyMOL: File > Run Script
# ============================================================

# Setup
bg_color white
set ray_opaque_background, on
hide everything
show cartoon, all

# Color all chains
color slate, chain R
color green, chain A
color lightpink, chain B
color paleyellow, chain C or chain G

# GPR55 surface
show surface, chain R
set surface_color, slate, chain R
set transparency, 0.7, chain R

# Create 3 HTVS target spheres
pseudoatom P0_ortho, pos=[110.02, 111.68, 81.10]
pseudoatom P3_allo, pos=[105.58, 102.05, 107.35]
pseudoatom Interface, pos=[112.678, 92.605, 112.545]

# Show spheres
show spheres, P0_ortho or P3_allo or Interface
set sphere_scale, 11, P0_ortho
set sphere_scale, 10, P3_allo
set sphere_scale, 12.5, Interface
set sphere_transparency, 0.3

# Colors
color cyan, P0_ortho
color tv_green, P3_allo
color hotpink, Interface

# Orient and render
orient
zoom all, 5
ray 1600, 1200
png D:/Antigravity AI/Major Project/Figures/PyMOL_Renders/HTVS_Target_Sites.png, dpi=600
