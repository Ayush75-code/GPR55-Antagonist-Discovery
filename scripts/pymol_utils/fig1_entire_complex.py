# Figure 1: Entire Protein Complex
# Run in PyMOL: File > Run Script

from pymol import cmd

OUTPUT = "D:/Antigravity AI/Major Project/Figures/PyMOL_Renders"

# Setup
cmd.bg_color("white")
cmd.hide("everything")
cmd.show("cartoon", "all")

# Consistent coloring
cmd.color("slate", "chain R")       # GPR55 - Blue
cmd.color("green", "chain A")       # Gα - Green
cmd.color("lightpink", "chain B")   # Gβ - Pink
cmd.color("paleyellow", "chain C or chain G")  # Gγ - Yellow

cmd.set("cartoon_transparency", 0.0)
cmd.set("ray_shadows", 0)
cmd.set("antialias", 2)

# Orient and render
cmd.orient()
cmd.zoom("all", 5)
cmd.ray(1600, 1200)
cmd.png(f"{OUTPUT}/Fig1_Entire_Complex.png", dpi=300)

print("✓ Saved: Fig1_Entire_Complex.png")
