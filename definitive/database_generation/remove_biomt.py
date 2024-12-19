#!/user/bin/env python

# Usage: python script.py path_to_PDBs_directory list

import os
import sys
import re
import pandas as pd
import subprocess
import traceback
from Bio.PDB import *
parser = PDBParser(QUIET = True)
import numpy
from math import sqrt
#from pymol import *


# ARGUMENTS:
pdb_list = sys.argv[2]
pdb_list = pd.read_csv(pdb_list, sep="\t", header=0)
fullPDB = sys.argv[1]
dir = os.getcwd()
selected = open("selected_pdb.csv", "w")
discarded = open("discarded_pdb.csv", "w")
errors = open("errors_pdb.csv", "w")

# EXECUTION:
counter = 0
for i in pdb_list.index:
    counter += 1
    pdb = pdb_list["pdb"][i]
    ent_file = os.path.join(fullPDB,f"pdb{pdb}.ent")
    try:
        pattern = r"REMARK 350   BIOMT1   2"
        with open(ent_file, "r") as file:
            file = file.read()
            match = re.search(pattern, file)
            print(match)
            if not match:
                selected.write(f"{pdb}\n")
            else:
                discarded.write(f"{pdb}\n")
            print(f"{pdb.upper()} has been classified correctly. Stucture number: {counter}")
    except:
        print(f"ERROR: {pdb.upper()} could not be classified. Stucture number: {counter}")
        errors.write(f"{pdb}\n")