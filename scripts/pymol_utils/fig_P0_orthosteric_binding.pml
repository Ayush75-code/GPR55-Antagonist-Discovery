# ============================================================
# Figure: AM251 Control at P0 Orthosteric Site (-9.458 kcal/mol)
# Run in PyMOL: File > Run Script
# ============================================================

# Load structures
load D:/Antigravity AI/Major Project/Phase 1 calidation (Control)/GPR55_receptor.pdbqt
load D:/Antigravity AI/Major Project/Phase 1 calidation (Control)/result/best_poses/overall_best/AM251_control_BEST_P0_-9.458kcalmol.pdbqt

# Setup
bg_color white
set ray_opaque_background, on
hide everything

# Show protein as cartoon
show cartoon, GPR55_receptor
color slate, GPR55_receptor

# Select ligand and pocket
select lig, AM251_control_BEST_P0_-9.458kcalmol
select pocket, byres GPR55_receptor within 4.5 of lig

# Pocket residues (Green Carbons)
show sticks, pocket
color green, pocket and elem C
color red, pocket and elem O
color blue, pocket and elem N
color yellow, pocket and elem S
set stick_radius, 0.15, pocket

# Ligand (Orange Carbons + spheres)
show sticks, lig
show sphere, lig
color orange, lig and elem C
color red, lig and elem O
color blue, lig and elem N
color yellow, lig and elem S
color green, lig and elem Cl
set sphere_scale, 0.25, lig
set stick_radius, 0.2, lig

# Add mesh around ligand
show mesh, lig
set mesh_color, orange, lig

# Residue labels
label pocket and name CA, "%s%s" % (resn, resi)
set label_color, black
set label_size, -0.5
set label_font_id, 7

# Final view and render
zoom (pocket or lig)
ray 1600, 1200
png D:/Antigravity AI/Major Project/Figures/PyMOL_Renders/AM251_Control_P0_Orthosteric.png, dpi=600
