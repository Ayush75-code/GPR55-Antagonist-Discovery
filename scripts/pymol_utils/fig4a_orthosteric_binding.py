# Figure 4a: Orthosteric Binding Site (P0) with Ligand
# Run in PyMOL: File > Run Script
# NOTE: Load your docked ligand PDBQT file first!

from pymol import cmd
from pymol.cgo import *
import math

OUTPUT = "D:/Antigravity AI/Major Project/Figures/PyMOL_Renders"

def make_mesh_sphere(name, center, radius, color, segments=30):
    x, y, z = center
    r, g, b = color
    obj = [LINEWIDTH, 4.0, COLOR, r, g, b]
    
    for i in range(segments + 1):
        theta = math.pi * i / segments
        obj.append(BEGIN)
        obj.append(LINE_STRIP)
        for j in range(segments + 1):
            phi = 2 * math.pi * j / segments
            px = x + radius * math.sin(theta) * math.cos(phi)
            py = y + radius * math.sin(theta) * math.sin(phi)
            pz = z + radius * math.cos(theta)
            obj.extend([VERTEX, px, py, pz])
        obj.append(END)
    
    for j in range(segments + 1):
        phi = 2 * math.pi * j / segments
        obj.append(BEGIN)
        obj.append(LINE_STRIP)
        for i in range(segments + 1):
            theta = math.pi * i / segments
            px = x + radius * math.sin(theta) * math.cos(phi)
            py = y + radius * math.sin(theta) * math.sin(phi)
            pz = z + radius * math.cos(theta)
            obj.extend([VERTEX, px, py, pz])
        obj.append(END)
    
    cmd.load_cgo(obj, name)

# Setup
cmd.bg_color("white")
cmd.hide("everything")
cmd.show("cartoon", "chain R")
cmd.color("slate", "chain R")
cmd.set("cartoon_transparency", 0.0)

# Pocket mesh (Orange)
make_mesh_sphere("pocket_P0", [110.02, 111.68, 81.10], 11, [1.0, 0.4, 0.0], 30)

# Show binding site residues
cmd.select("site_P0", "chain R within 5 of (x=110 and y=111 and z=81)")
cmd.show("sticks", "site_P0")

# Show ligand as ball-and-stick with mesh
cmd.show("sticks", "organic")
cmd.show("spheres", "organic")
cmd.set("sphere_scale", 0.25, "organic")
cmd.set("stick_radius", 0.12, "organic")

# Ligand coloring by element
cmd.color("cyan", "organic and elem C")
cmd.color("red", "organic and elem O")
cmd.color("blue", "organic and elem N")
cmd.color("yellow", "organic and elem S")
cmd.color("green", "organic and elem Cl")
cmd.color("orange", "organic and elem Br")
cmd.color("white", "organic and elem H")

# Render
cmd.set("ray_shadows", 0)
cmd.center("site_P0")
cmd.zoom("site_P0", 8)
cmd.ray(1600, 1200)
cmd.png(f"{OUTPUT}/Fig4a_Orthosteric_P0_Binding.png", dpi=300)

print("✓ Saved: Fig4a_Orthosteric_P0_Binding.png")
