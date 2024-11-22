#!/user/bin/env python

# Usage: python script.py fullPDB_path PDB_list

import os
import sys
import pandas as pd
import subprocess
import traceback
from Bio.PDB import *
parser = PDBParser(QUIET=True)

# LISTS:
#oligomerization_states = ["monomeric", "dimeric", "trimeric", "tetrameric", "pentameric", "hexameric",
                          #"heptameric", "octameric","nonameric"]
oligomerization_states = {
	"monomeric": 1,
	"dimeric": 2,
    "trimeric":3,
    "tetrameric":4,
    "pentameric":5,
    "hexameric":6,
    "heptameric":7,
    "octameric":8,
    "nonameric":9
}

# ARGUMENTS:
pdb_list = sys.argv[2]
pdb_list = pd.read_csv(pdb_list, sep="\t", header=0)
#pdb_list = pdb_list.set_index(pdb_list[0])
fullPDB = sys.argv[1]
# annotated = open("annotated.csv","w")
# annotated.write(f"pdb\tannotation\n")
non_annotated = open("non_annotated.csv", "w")
non_annotated.write(f"pdb\n")

# EXECUTION:
for i in pdb_list.index:
    pdb = pdb_list.loc[i,"pdb"]
    ent_file = os.path.join(fullPDB,f"pdb{pdb}.ent")
    try:
        structure = parser.get_structure(pdb, ent_file)
        chains = structure.get_chains()
        chain_number = 0
        for chain in chains:
            chain_number += 1
        pdb_list.at[i, "chain_number"] = chain_number
    except:
        print(f"Error in {pdb}: impossible to annotate chain number.")
        non_annotated.write(f"{pdb}\n")
        #pdb_list.at[i, "chain_number"] = ""
    try:
        result = subprocess.check_output(f"grep 'AUTHOR DETERMINED BIOLOGICAL UNIT' {ent_file}", shell=True)
        result = str(result)
        for state in oligomerization_states.keys():
            if state.upper() in result:
                oligomer = oligomerization_states[state]
                pdb_list.at[i, "oligomer"] = oligomer
    except:
        non_annotated.write(f"{pdb}\n")
        oligomer=""
        #pdb_list.at[i, "oligomer"] = ""

non_annotated.close()
pdb_list.to_csv("annotated_chains_oligomer.csv",sep="\t")

'''
CHECK FOR CASES WITH MORE THAN ONE BIOMOLECULE WITH IRA
'''