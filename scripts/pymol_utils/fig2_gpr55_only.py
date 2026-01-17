# Figure 2: GPR55 Receptor Only
# Run in PyMOL: File > Run Script

from pymol import cmd

OUTPUT = "D:/Antigravity AI/Major Project/Figures/PyMOL_Renders"

# Setup
cmd.bg_color("white")
cmd.hide("everything")
cmd.show("cartoon", "chain R")  # Only GPR55

# Color GPR55
cmd.color("slate", "chain R")

cmd.set("cartoon_transparency", 0.0)
cmd.set("ray_shadows", 0)
cmd.set("antialias", 2)

# Orient and render
cmd.orient("chain R")
cmd.zoom("chain R", 5)
cmd.ray(1600, 1200)
cmd.png(f"{OUTPUT}/Fig2_GPR55_Only.png", dpi=300)

print("✓ Saved: Fig2_GPR55_Only.png")
