#!/user/bin/env python

# Usage: python script.py path_to_PDB list errors_list_name

import os
import sys
import traceback
import subprocess

# ARGUMENTS:
fullPDB = sys.argv[1]
pdb_list = sys.argv[2]
pdb_list = open(pdb_list, mode='r').read().splitlines()
errors_name = sys.argv[3]
errors = open(f"{errors_name}.txt", "w")
dir = os.getcwd()

# EXECUTION:
if "output_annotation" not in os.listdir(dir):
    subprocess.run(f"mkdir output_annotation", shell=True)
counter = 0
errors_num = 0
for pdb in pdb_list:
    counter+=1
    try:
        ent_file = os.path.join(fullPDB,f"pdb{pdb}.ent")
        if f"{pdb}" not in os.listdir(f"{dir}/output_annotation"):
            subprocess.run(f"mkdir output_annotation/{pdb}", shell=True)
        if f"pdb{pdb}_input_search.pdb" not in os.listdir(f"{dir}/output_annotation/{pdb}/"):
            subprocess.run(f"python /cri4/iracema/REPO_GIT/borges-arcimboldo/ARCIMBOLDO_FULL/ALEPH/aleph/core/ALEPH.py"
                           f" annotate --pdbmodel {ent_file} --strictness_ah 0.5 --strictness_bs 0.3", shell=True)
            subprocess.run(f"mv pdb* strictnesses.pdb output.json fromCVtoAA_12.txt CVs.txt output_annotation/{pdb}", shell=True)
        print(f"{pdb.upper()} has been correctly annotated. Structure number {counter}.")
    except:
        print(f"Error in {pdb.upper()}. Could not be annotated. Structure number {counter}.")
        traceback.print_exc()
        errors_num += 1
        errors.write(f"{pdb}\n")

print(f"\n")
print(f"OUTPUT SUMMARY:")
print(f"Total correctly annotated structures: {counter-errors_num}.")
print(f"Total failed annotations: {errors_num}.")
