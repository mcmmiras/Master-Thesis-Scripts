#!/usr/bin/env python

# Usage: open pdb(file) in pymol and manually execute the script

# Importation of the library:
from pymol import *

# Setting the spectrum parameters:
#cmd.spectrum(expression = "b",palette = "rainbow_rev",selection = "all",minimum = 0,maximum = 100) # for normal settings
cmd.spectrum(expression = "b",palette = "blue_red",selection = "all",minimum = 0,maximum = 100) # for heptads detection
# Setting the cartoon oval length:
cmd.set(name="cartoon_oval_length", value=1)

# Setting the background color:
cmd.bg_color("white")

# Setting the shadows:
cmd.set('ray_shadows','off')