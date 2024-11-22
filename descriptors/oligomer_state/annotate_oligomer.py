#!/user/bin/env python

# Usage: python script.py fullPDB_path PDB_list

import os
import sys
import pandas as pd
import subprocess
import traceback

# LISTS:
oligomerization_states = ["monomeric", "dimeric", "trimeric", "tetrameric", "pentameric", "hexameric",
                          "heptameric", "octameric","nonameric"]

# ARGUMENTS:
pdb_list = sys.argv[2]
pdb_list = pd.read_csv(pdb_list, sep="\t", header=0)
#pdb_list = pdb_list.set_index(pdb_list[0])
fullPDB = sys.argv[1]
annotated = open("annotated.csv","w")
annotated.write(f"pdb\tannotation\n")
non_annotated = open("non_annotated.csv", "w")
non_annotated.write(f"pdb\n")

# EXECUTION:
for i in pdb_list.index:
    pdb = pdb_list.loc[i,"pdb"]
    ent_file = os.path.join(fullPDB,f"pdb{pdb}.ent")
    try:
        result = subprocess.check_output(f"grep 'AUTHOR DETERMINED BIOLOGICAL UNIT' {ent_file}", shell=True)
        result = str(result)
        for state in oligomerization_states:
            state = state.upper()
            if state in result:
                oligomer = state
                annotated.write(f"{pdb}\t{oligomer}\n")
    except:
        non_annotated.write(f"{pdb}\n")
        oligomer=""
        annotated.write(f"{pdb}\t{oligomer}\n")

non_annotated.close()
annotated.close()

'''
CHECK FOR CASES WITH MORE THAN ONE BIOMOLECULE WITH IRA
'''