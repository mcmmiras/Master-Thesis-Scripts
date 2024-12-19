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
result_cc = open(f"confirmed_chain_helix.csv","w")
result_cc.write(f"pdb\tah_count\tah_res\ttotal_res\tah_res/total_res\t"
                f"chains\tdetermination\toligomer\toligomer_str\n")
result_ccOlig = open(f"confirmed_chain_helix_oligomer.csv","w")
result_ccOlig.write(f"pdb\tah_count\tah_res\ttotal_res\tah_res/total_res\t"
                f"chains\tdetermination\toligomer\toligomer_str\n")
result_notcc = open(f"nonconfirmed_oligomer.csv","w")
result_notcc.write(f"pdb\tah_count\tah_res\ttotal_res\tah_res/total_res\t"
                f"chains\tdetermination\toligomer\toligomer_str\n")
#pdb_list = pdb_list.set_index(pdb_list[0])
fullPDB = sys.argv[1]
# annotated = open("annotated.csv","w")
# annotated.write(f"pdb\tannotation\n")
non_annotated = open("non_annotated.csv", "w")
non_annotated.write(f"pdb\n")

# EXECUTION:
if f"annotated_chains_oligomer.csv" not in os.listdir(os.getcwd()):
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
                    pdb_list.at[i, "determination"] = "author"
                    pdb_list.at[i, "oligomer"] = oligomer
                    pdb_list.at[i, "oligomer_str"] = state
        except:
            try:
                result = subprocess.check_output(f"grep 'SOFTWARE DETERMINED QUATERNARY STRUCTURE:' {ent_file}", shell=True)
                result = str(result)
                for state in oligomerization_states.keys():
                    if state.upper() in result:
                        oligomer = oligomerization_states[state]
                        pdb_list.at[i,"determination"] = "software"
                        pdb_list.at[i, "oligomer"] = oligomer
                        pdb_list.at[i, "oligomer_str"] = state
            except:
                non_annotated.write(f"{pdb}\n")
                oligomer=""
                traceback.print_exc()
                #pdb_list.at[i, "oligomer"] = ""

    non_annotated.close()
    pdb_list.to_csv("annotated_chains_oligomer.csv",sep="\t")


result = pd.read_csv("annotated_chains_oligomer.csv", sep=f"\t")
for i in result.index:
    if result.loc[i,"helix_count"] == result.loc[i,"chain_number"]:
        result_cc.write(f"{result.loc[i,'pdb']}\t"
                        f"{result.loc[i,'helix_count']}\t"
                        f"{result.loc[i,'helix_residues']}\t"
                        f"{result.loc[i,'total_residues']}\t"
                        f"{result.loc[i,'helix_residues/total_residues']}\t"
                        f"{result.loc[i, 'chain_number']}\t"
                        f"{result.loc[i, 'determination']}\t"
                        f"{result.loc[i, 'oligomer']}\t"
                        f"{result.loc[i, 'oligomer_str']}\n")
        if result.loc[i,"chain_number"] == result.loc[i,"oligomer"]:
            result_ccOlig.write(f"{result.loc[i, 'pdb']}\t"
                            f"{result.loc[i, 'helix_count']}\t"
                            f"{result.loc[i, 'helix_residues']}\t"
                            f"{result.loc[i, 'total_residues']}\t"
                            f"{result.loc[i, 'helix_residues/total_residues']}\t"
                            f"{result.loc[i, 'chain_number']}\t"
                            f"{result.loc[i, 'determination']}\t"
                            f"{result.loc[i, 'oligomer']}\t"
                            f"{result.loc[i, 'oligomer_str']}\n")
    else:
        result_notcc.write(f"{result.loc[i, 'pdb']}\t"
                            f"{result.loc[i, 'helix_count']}\t"
                            f"{result.loc[i, 'helix_residues']}\t"
                            f"{result.loc[i, 'total_residues']}\t"
                            f"{result.loc[i, 'helix_residues/total_residues']}\t"
                            f"{result.loc[i, 'chain_number']}\t"
                            f"{result.loc[i, 'determination']}\t"
                            f"{result.loc[i, 'oligomer']}\t"
                            f"{result.loc[i, 'oligomer_str']}\n")

'''
CHECK FOR CASES WITH MORE THAN ONE BIOMOLECULE WITH IRA
'''