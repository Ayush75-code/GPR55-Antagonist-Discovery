# ============================================================
# GPR55 Project - Complete Visualization Script
# Generates all figures for presentation
# ============================================================
# Run in PyMOL: File > Run Script > select this file
# OR: @D:/Antigravity AI/Major Project/Scripts/PyMOL_Utils/generate_all_figures.py
# ============================================================

from pymol import cmd
from pymol.cgo import *
import math

# Output directory
OUTPUT_DIR = "D:/Antigravity AI/Major Project/Figures/PyMOL_Renders"

# ============================================================
# HELPER FUNCTION: Create Mesh Sphere
# ============================================================
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

# ============================================================
# COMMON SETUP
# ============================================================
def setup_display():
    cmd.bg_color("white")
    cmd.hide("everything")
    cmd.show("cartoon", "all")
    
    # Consistent coloring across all images
    cmd.color("slate", "chain R")       # GPR55 - Blue
    cmd.color("green", "chain A")       # Gα - Green
    cmd.color("lightpink", "chain B")   # Gβ - Pink
    cmd.color("paleyellow", "chain C or chain G")  # Gγ - Yellow
    
    cmd.set("cartoon_transparency", 0.0)
    cmd.set("ray_shadows", 0)
    cmd.set("antialias", 2)

def render_image(filename, width=1600, height=1200, dpi=300):
    cmd.ray(width, height)
    cmd.png(f"{OUTPUT_DIR}/{filename}", dpi=dpi)
    print(f"✓ Saved: {filename}")

# ============================================================
# FIGURE 1: Entire Protein Complex
# ============================================================
def figure1_entire_complex():
    print("\n=== Figure 1: Entire Protein Complex ===")
    setup_display()
    cmd.orient()
    cmd.zoom("all", 5)
    render_image("Fig1_Entire_Complex.png")

# ============================================================
# FIGURE 2: GPR55 Only
# ============================================================
def figure2_gpr55_only():
    print("\n=== Figure 2: GPR55 Receptor Only ===")
    setup_display()
    cmd.hide("cartoon", "chain A or chain B or chain C or chain G")
    cmd.orient("chain R")
    cmd.zoom("chain R", 5)
    render_image("Fig2_GPR55_Only.png")
    # Restore
    cmd.show("cartoon", "all")

# ============================================================
# FIGURE 3: All 6 Binding Pockets (Phase 1)
# ============================================================
def figure3_six_pockets():
    print("\n=== Figure 3: Six Binding Pockets ===")
    setup_display()
    
    # Define all 6 pockets
    pockets = {
        "P0": {"center": [110.02, 111.68, 81.10], "radius": 11, "color": [1.0, 0.5, 0.0]},    # Orange
        "P1": {"center": [122.01, 102.28, 106.69], "radius": 10, "color": [0.8, 0.8, 0.0]},   # Yellow
        "P2": {"center": [97.23, 106.29, 104.49], "radius": 10, "color": [0.6, 0.0, 0.6]},    # Purple
        "P3": {"center": [105.58, 102.05, 107.35], "radius": 10, "color": [0.0, 0.6, 0.0]},   # Green
        "P4": {"center": [111.03, 95.74, 101.59], "radius": 10, "color": [0.0, 0.6, 0.8]},    # Cyan
        "P5": {"center": [111.95, 114.75, 103.41], "radius": 10, "color": [0.8, 0.4, 0.6]},   # Pink
    }
    
    # Create mesh spheres for each pocket
    for name, data in pockets.items():
        make_mesh_sphere(f"pocket_{name}", data["center"], data["radius"], data["color"], 25)
        
        # Add center point
        cmd.pseudoatom(f"center_{name}", pos=data["center"])
        cmd.show("spheres", f"center_{name}")
        cmd.set("sphere_scale", 1.0, f"center_{name}")
    
    # Color center points
    cmd.color("orange", "center_P0")
    cmd.color("yellow", "center_P1")
    cmd.color("purple", "center_P2")
    cmd.color("green", "center_P3")
    cmd.color("cyan", "center_P4")
    cmd.color("pink", "center_P5")
    
    cmd.orient("chain R")
    cmd.zoom("chain R", 8)
    render_image("Fig3_Six_Binding_Pockets.png")
    
    # Cleanup
    for name in pockets.keys():
        cmd.delete(f"pocket_{name}")
        cmd.delete(f"center_{name}")

# ============================================================
# FIGURE 4a: Orthosteric Site (P0) with Ligand
# ============================================================
def figure4a_orthosteric_binding():
    print("\n=== Figure 4a: Orthosteric Binding (P0) ===")
    setup_display()
    
    # Show only GPR55
    cmd.hide("cartoon", "chain A or chain B or chain C or chain G")
    
    # Create pocket mesh
    make_mesh_sphere("pocket_P0", [110.02, 111.68, 81.10], 11, [1.0, 0.4, 0.0], 30)
    
    # If ligand is loaded, show it as ball-and-stick with mesh
    # Assuming ligand is named "AM251" or similar
    try:
        cmd.show("sticks", "organic")
        cmd.show("spheres", "organic")
        cmd.set("sphere_scale", 0.3, "organic")
        cmd.set("stick_radius", 0.15, "organic")
        cmd.color("cyan", "organic and elem C")
        cmd.color("red", "organic and elem O")
        cmd.color("blue", "organic and elem N")
        cmd.color("yellow", "organic and elem S")
        cmd.color("green", "organic and elem Cl")
    except:
        pass
    
    # Show binding site residues as sticks
    cmd.select("site_P0", "chain R within 5 of (x=110 and y=111 and z=81)")
    cmd.show("sticks", "site_P0")
    
    # Zoom to orthosteric site
    cmd.center("site_P0")
    cmd.zoom("site_P0", 10)
    
    render_image("Fig4a_Orthosteric_P0_Binding.png")
    
    # Cleanup
    cmd.delete("pocket_P0")
    cmd.delete("site_P0")
    cmd.show("cartoon", "all")

# ============================================================
# FIGURE 4b: Allosteric Site (P3) with Ligand
# ============================================================
def figure4b_allosteric_binding():
    print("\n=== Figure 4b: Allosteric Binding (P3) ===")
    setup_display()
    
    # Show only GPR55
    cmd.hide("cartoon", "chain A or chain B or chain C or chain G")
    
    # Create pocket mesh
    make_mesh_sphere("pocket_P3", [105.58, 102.05, 107.35], 10, [0.0, 0.6, 0.0], 30)
    
    # Show ligand if present
    try:
        cmd.show("sticks", "organic")
        cmd.show("spheres", "organic")
        cmd.set("sphere_scale", 0.3, "organic")
        cmd.set("stick_radius", 0.15, "organic")
        cmd.color("magenta", "organic and elem C")
    except:
        pass
    
    # Show binding site residues as sticks
    cmd.select("site_P3", "chain R within 5 of (x=105 and y=102 and z=107)")
    cmd.show("sticks", "site_P3")
    
    # Zoom to allosteric site
    cmd.center("site_P3")
    cmd.zoom("site_P3", 10)
    
    render_image("Fig4b_Allosteric_P3_Binding.png")
    
    # Cleanup
    cmd.delete("pocket_P3")
    cmd.delete("site_P3")
    cmd.show("cartoon", "all")

# ============================================================
# MAIN: Generate All Figures
# ============================================================
def generate_all():
    print("\n" + "="*60)
    print("  GPR55 Project - Generating All Presentation Figures")
    print("="*60)
    
    figure1_entire_complex()
    figure2_gpr55_only()
    figure3_six_pockets()
    figure4a_orthosteric_binding()
    figure4b_allosteric_binding()
    
    print("\n" + "="*60)
    print(f"  All figures saved to: {OUTPUT_DIR}")
    print("="*60 + "\n")

# Run if executed as script
if __name__ == "__main__":
    generate_all()
else:
    print("\nScript loaded! Run: generate_all()")
    print("Or run individual figures:")
    print("  figure1_entire_complex()")
    print("  figure2_gpr55_only()")
    print("  figure3_six_pockets()")
    print("  figure4a_orthosteric_binding()")
    print("  figure4b_allosteric_binding()")
