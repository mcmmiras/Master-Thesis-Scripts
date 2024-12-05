#!/usr/bin/env python

# Usage: pymol script.py folder

# Importation of the library:
from pymol import *
import sys, os

counter=0
# Loading the structure:
for ele in os.listdir(sys.argv[1]):
    pdb = ele.split("_")[1].split("-")[0]
    scale = (ele.split("-")[1]).split(".")[0]
    cmd.load(os.path.join(sys.argv[1],ele))
    # Setting the spectrum parameters:
    cmd.spectrum(expression = "b",palette = "blue_green",selection = "all",minimum = 0,maximum = 100)
    # Setting the cartoon oval length:
    cmd.set(name="cartoon_oval_length", value=1)
    # Setting the background color:
    cmd.bg_color("white")
    # Setting the shadows:
    cmd.set('ray_shadows','off')
    # Saving the image:
    cmd.draw()
    cmd.save(f"colored_hydrophobicity_{pdb}_{scale}.png")
    counter += 1
    print(f"{pdb.upper()} structure with applied {scale} hydrophobicity scale has been saved.")
    cmd.delete(all)
print(f"All structures finished. Total colored structures: {counter}.")
cmd.quit()