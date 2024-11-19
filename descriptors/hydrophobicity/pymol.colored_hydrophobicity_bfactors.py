#!/usr/bin/env python

# Usage: pymol script.py pbd/s(file)

# Importation of the library:
from pymol import *
# Loading the structure:
#for i in range(2,len(sys.argv)):
    # Loading the PDB file:
    #cmd.load(f"{sys.argv[i]}")
# Setting the spectrum parameters:
cmd.spectrum(expression = "b",palette = "red_white_blue",selection = "all",minimum = 0,maximum = 100)
# Setting the cartoon oval length:
cmd.set(name="cartoon_oval_length", value=1)
# Setting the background color:
cmd.bg_color("white")
# Setting the shadows:
cmd.set('ray_shadows','off')