#!/user/bin/env python

# Usage: python script.py path_to_fullPDB pdb_list(.txt)
# This script detects which helical structures are considered coiled-coils.


# INFORMATION:

# SOCKET:
# What you need to run SOCKET:
# 	1)	A PDB-format file to examine (ex., ent file)
# 	2)	A DSSP-format file generated from the PDB file

# DSSP:
# IMPORTANT: If you are using the recent version of DSSP, you must cut the 1-136 (both included) of DSSP output file
# for further use with SOCKET 2.0.
#   1)  This is a rewrite of DSSP, now offering full mmCIF support. The difference with previous releases of DSSP
#       is that it now writes out an annotated mmCIF file by default, storing the secondary structure information
#       in the _struct_conf category.
#   2)  Another new feature in this version of DSSP is that it now defines Poly-Proline helices as well.


# LIBRARY IMPORTATION:
import subprocess
import sys
import pandas as pd
import os


# ARGUMENTS:
counter = 0
fullPDB = sys.argv[1]
pdb_list = sys.argv[2]
pdb_list = open(pdb_list, "r")
pdb_list = pd.read_csv(pdb_list,sep='\t', usecols=[0], header = None)
pdb_list = pdb_list[0].tolist()


# OUTPUT FILES:
confirmedCC = open("socket_validated_coiled_coils.txt", "w")
notCC = open("socket_validated_NOT_coiled_coils.txt", "w")
unclassified = open("socket_unclassified.txt","w")
errorPDB = open("socket_ERROR.txt", "w")

# DIRECTORIES CREATION:
dir = os.getcwd()
if not "output_dssp_files" in os.listdir(dir):
    subprocess.run(f"mkdir output_dssp_files", shell=True)
if not "output_socket_files" in os.listdir(dir):
    subprocess.run(f"mkdir output_socket_files", shell=True)


# FUNCTIONS:
# Analysis function: starting from an ent file, it creates an output after DSSP and SOCKET2 runs.
def analysis(ent, dssp, out): # Faltaría final_dssp,
    subprocess.run("mkdssp -v " + ent + " " + dssp, shell=True)  # Create dssp file
    #subprocess.run("cut -c 1-136 " + dssp + " > " + final_dssp, shell=True) # Cut chrs 1-136 from dssp file
    subprocess.run("/cri4/mari/repositories/socket2_linux -f " + ent + " -s " + dssp + ">" + out, shell=True) # Run Socket2
    subprocess.run("grep result " + out, shell=True)

# EXECUTION:
# Classification of PDBs from the provided txt file as coiled-coils or NOT coiled-coils:
for pdb in pdb_list:
    try:
        counter += 1
        ent_file = f"{fullPDB}{pdb}.ent"
        dssp_file = f"./output_dssp_files/{pdb}.dssp"
        #final_dssp_file = "cut_" + dssp_file
        out_file = f"./output_socket_files/result_{pdb}.out"
        try:
            analysis(ent_file, dssp_file, out_file)
            with open(out_file,"r") as result:
                result = result.read()
                if "0 result NO COILED COILS" in result:
                    notCC.write(f"{pdb}\n")
                elif "COILED COILS PRESENT" in result:
                    confirmedCC.write(f"{pdb}\n")
                else:
                    unclassified.write(f"{pdb}\n")
            print(f"{pdb} has been classified. Structure number: {counter}.\n\n")
        except:
            print(f"Error in {pdb}: SOCKET cannot compute analysis.")
            errorPDB.write(f"{pdb}\n")
    except:
        print(f"Error in {pdb}: not found.")

confirmedCC.close()
notCC.close()
unclassified.close()
errorPDB.close()