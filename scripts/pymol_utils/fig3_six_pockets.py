# Figure 3: Six Binding Pockets (Phase 1 Control)
# Run in PyMOL: File > Run Script

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
cmd.show("cartoon", "chain R")  # Only GPR55
cmd.color("slate", "chain R")
cmd.set("cartoon_transparency", 0.0)

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
    cmd.pseudoatom(f"center_{name}", pos=data["center"])
    cmd.show("spheres", f"center_{name}")
    cmd.set("sphere_scale", 1.0, f"center_{name}")

# Color centers
cmd.color("orange", "center_P0")
cmd.color("yellow", "center_P1")
cmd.color("purple", "center_P2")
cmd.color("green", "center_P3")
cmd.color("cyan", "center_P4")
cmd.color("pink", "center_P5")

# Render
cmd.set("ray_shadows", 0)
cmd.orient("chain R")
cmd.zoom("chain R", 8)
cmd.ray(1600, 1200)
cmd.png(f"{OUTPUT}/Fig3_Six_Binding_Pockets.png", dpi=300)

print("✓ Saved: Fig3_Six_Binding_Pockets.png")
