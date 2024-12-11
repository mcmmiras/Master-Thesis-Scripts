#!/user/bin/env python

# Usage: python script.py path_to_PDB list

import os
import sys
import pandas as pd
import subprocess

# ARGUMENTS:
fullPDB = sys.argv[1]
pdb_list = sys.argv[2]
pdb_list = open(pdb_list, mode='r').read().splitlines()
dir = os.getcwd()
errors = open("errors_pdb.csv", "w")

# EXECUTION:
subprocess.run(f"mkdir output_annotation", shell=True)
counter = 0
list_disc_res = list()
for pdb in pdb_list:
    counter+=1
    ent_file = os.path.join(fullPDB,f"pdb{pdb}.ent")
    subprocess.run(f"mkdir output_annotation/{pdb}", shell=True)
    subprocess.run(f"python /cri4/iracema/REPO_GIT/borges-arcimboldo/ARCIMBOLDO_FULL/ALEPH/aleph/core/ALEPH.py annotate "
                   f"--pdbmodel {ent_file} --strictness_ah 0.5 --strictness_bs 0.3", shell=True)
    subprocess.run(f"mv pdb* strictnesses.pdb output.json fromCVtoAA_12.txt CVs.txt output_annotation/{pdb}", shell=True)
    print(f"{pdb.upper()} has been correctly annotated. Structure number {counter}.")
    #except:
    #print(f"Error in {pdb.upper()}. Could not be annotated.")
    #errors.write(f"{pdb}\n")
