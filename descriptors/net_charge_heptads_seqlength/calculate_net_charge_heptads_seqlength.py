# !usr/bin/env python
from enum import unique

# Usage: python script.py fullPDB_path PDB_list(.txt)


# LIBRARIES IMPORTATION:
from Bio.SeqUtils import IsoelectricPoint
from Bio.PDB import *
parser = PDBParser(QUIET=True)
io = PDBIO()
sr = ShrakeRupley()
import os, sys, subprocess
import pandas as pd
import re
import traceback


# ARGUMENTS:
fullPDB = sys.argv[1]
list = sys.argv[2]
df = pd.read_csv(list, sep="\t")
dir = os.getcwd()
if not "fasta_files" in os.listdir(dir):
    subprocess.run(f"mkdir fasta_files/", shell=True)
if not "modified_PDBs" in os.listdir(dir):
    subprocess.run(f"mkdir modified_PDBs/", shell=True)
if not "heptads_annotation" in os.listdir(dir):
    subprocess.run(f"mkdir heptads_annotation/", shell=True)
resultsHeptads = open("output_heptads_annotation.txt", "w")
errors = open("unprocessed_files_traceback.txt","w")


# DICTIONARIES:
helicalStructure = {
#    "hydrophobic": ["A", "C", "G", "I", "L","M","F","P","W","Y","V"],
#    "hydrophilic": ["N", "Q", "S", "T"],
#    "acid": ["D", "E"],
#    "base": ["R", "H", "K"]
#    "hydrophobic": ["A", "I", "L", "M", "F", "V", "P", "G"],
    "hydrophobic": ["A", "I", "L", "V"], # Most common hydrophobic residues found in a/d cores from helices.
    "hydrophilic": ["N", "Q", "S", "T"],
    "charged": ["R", "K", "D", "E"],
    "polar": ["Q", "N", "H", "S", "T", "Y", "C"],
    "amphipathic": ["W", "Y", "M"]
}

aminoacids = {
    "ALA": "A",
    "CYS": "C",
    "ASP": "D",
    "GLU": "E",
    "PHE": "F",
    "GLY": "G",
    "HIS": "H",
    "ILE": "I",
    "LYS": "K",
    "LEU": "L",
    "MET": "M",
    "ASN": "N",
    "PRO": "P",
    "GLN": "Q",
    "ARG": "R",
    "SER": "S",
    "THR": "T",
    "VAL": "V",
    "TRP": "W",
    "TYR": "Y"
}


# FUNCTIONS:
# 1. Function to assign heptad repeats:
def assign_heptads(pdb, index):
    df_heptads = pd.DataFrame()
    heptads = str()
    #hydrophobicResidues = str()
    if not f"{pdb}.fasta" in os.listdir(f"{dir}/fasta_files/"):
        subprocess.run(f"wget -O ./fasta_files/{pdb}.fasta https://www.rcsb.org/fasta/entry/{pdb}", shell=True)
    fasta_file = f"{dir}/fasta_files/{pdb}.fasta"
    fasta_file = open(fasta_file, "r")
    fasta = fasta_file.read().splitlines()
    df_heptads.at[index, f"pdb"] = pdb
    unique_chains = 0
    resultsHeptads.write(f"{pdb}\n")
    for line in fasta:
        if ">" not in line:
            unique_chains += 1
    df_heptads.at[index,"unique_chains"] = unique_chains
    i_seq = 0
    for line in fasta:
        # for i_seq, seq in enumerate(sequences_list):
        #     seq_num = i_seq + 1
        #     df_heptads[f"seq1_heptad{heptad_num}_position"] = pos
        i = 0
        if ">" not in line:
            i_seq += 1
            seq = line
            seq_length = int(len(seq))
            while i < (len(seq)):
                if seq[i] == "X":
                    i+=1
                    heptads = heptads+"-"
                else:
                    if (seq[i] in helicalStructure["hydrophobic"] and seq[i+3] in helicalStructure["hydrophobic"]):
                        heptads = heptads+"abcdefg"
                        i+=7
                    else:
                        i+=1
                        heptads = heptads+"-"
            heptads = heptads[:len(seq)]
            positions = [i.start() for i in re.finditer("abcdefg", heptads)]
            resultsHeptads.write(f"{seq}\n{heptads}\n{positions}\n\n")
            df_heptads.at[index, f"sequence{str(i_seq)}_length"] = seq_length
            df_heptads.at[index, f"seq{str(i_seq)}_heptad_number"] = len(positions)
            for i_pos, pos in enumerate(positions):
                heptad_num = i_pos+1
                df_heptads.at[index,f"seq{str(i_seq)}_heptad{str(heptad_num)}_position"] = pos
        # for ele in heptads:
        #     if ele == "a" or ele == "d":
        #         pos = heptads.index(ele)
        #         hydrophobicResidues = hydrophobicResidues+seq[pos]
            #return seq, positions, seq_length, unique_chains#, hydrophobicResidues
    fasta_file.close()
    return df_heptads


# 2. Function to calculate the net-charge of a protein in a specific pH.
#   - Used pH: 7 (physiological)
def calculate_net_charge(ent):
    IP = IsoelectricPoint.IsoelectricPoint(ent)
    charge = round(IP.charge_at_pH(7.0), 2)
    return charge


# 3. Function to change b-factors following heptad annotation:
# Saving results: Save the modified structure to a new PDB file:
def bfactorToHeptadAnnotation(structure, pdb, pos, seq):
    residuesList = []
    residuesElements = []
    abvResiduesList = str()
    foundHeptads = []
    for model in structure:
        for chain in model:
            if chain.id == "A":
                residues = chain.get_residues()
                for residue in residues:
                    residuesElements.append(residue)
                    if residue.get_resname() != "HOH":
                        residuesList.append(residue.get_resname())
    for position in pos:
        foundHeptads.append(seq[position:position+7])
    for residue in residuesList:
        if residue in aminoacids:
            abvResiduesList = abvResiduesList + aminoacids[residue]
        else:
            abvResiduesList = abvResiduesList + "-"
    residues_set = set()  # Track residues with modified B-factors
    for heptad in foundHeptads:
        start_index = 0
        while (location := abvResiduesList.find(heptad, start_index)) != -1:
            for i in range(location, location + 7):
                residues_set.add(residuesElements[i])
                for atom in residuesElements[i]:
                    atom.set_bfactor(100 if i % 7 == 0 or i % 7 == 4 else 50)
            start_index = location + 1
    for residue in residuesElements:
        if residue not in residues_set:
            for atom in residue:
                atom.set_bfactor(0)
    io.set_structure(structure)
    #if not f"modified_{pdb}.pdb" in os.listdir(f"{dir}/modified_PDBs/"):
    io.save(f"{dir}/modified_PDBs/modified_{pdb}.pdb")
    return foundHeptads



# EXECUTION:
# Iteration through all the PDBs list:
for i in df.index:
    try: # Obtain PDB from list
        pdb = df['pdb'][i]
        ent_file = f"pdb{df['pdb'][i]}.ent"
        ent_file = os.path.join(fullPDB, ent_file)
        try: # Calculate protein's global net-charge
            try: # Assign heptad repeats
                structure = parser.get_structure(pdb,ent_file)
                #sequence, positions, seq_length, unique_chains =
                df_heptads = assign_heptads(pdb,i)
                df_heptads.at[i,"net_charge"] = calculate_net_charge(ent_file)
                df_heptads.to_csv(f"heptads_annotation/{pdb}_heptads_annotation.csv", sep="\t")
                #df.at[i, "sequence"] = sequence
                #df.at[i, "unique_chains"] = unique_chains
                #df.at[i, "sequence_length"] = seq_length
                #if len(positions) > 0:
                 #   df.at[i, "heptads"] = "yes"
                #else:
                 #   df.at[i, "heptads"] = "no"
                #df.at[i, "heptad_number"] = len(positions)
                #df.at[i, "heptad_positions"] = str(positions)
                ########################3# sr.compute(structure, level="S")
                # print(round(structure.sasa,2))
                # if hydrophobicCore:
                #     try: # Calculate hydrophobic core's net-charge
                #         df.at[i, "hydrophobic_core(HC)"] = hydrophobicCore
                #         df.at[i, "HC_net_charge"] = calculate_net_charge(hydrophobicCore)
                #     except:
                #         print(f"Error in {pdb}: impossible to calculate the net charge.")
                # try:
                #     heptadRepeats = bfactorToHeptadAnnotation(structure, pdb, positions, sequence)
                #     df.at[i, "heptad_repeats"] = ",".join(heptadRepeats)
                # except:
                #     print(f"Error in {pdb}: impossible to modify b-factors.")
                #     errors.write(f"{pdb}\nError in {pdb}: impossible to modify b-factors.\n")
                #     traceback.print_exc(file=errors)
                #     errors.write(f"\n\n")
            except:
                print(f"Error in {pdb}: impossible to assign heptad repeats.")
                errors.write(f"{pdb}\nError in {pdb}: impossible to assign heptad repeats.\n")
                traceback.print_exc(file=errors)
                errors.write(f"\n\n")
        except:
            print(f"Error in {pdb}: impossible to calculate global net charge.")
            errors.write(f"{pdb}\nError in {pdb}: impossible to calculate global net charge.\n")
            traceback.print_exc(file=errors)
            errors.write(f"\n\n")
    except:
        print(f"Error: {pdb} PDB not found!")
        errors.write(f"{pdb}\nError: {pdb} PDB not found!\n")
        traceback.print_exc(file=errors)
        errors.write(f"\n\n")

resultsHeptads.close()
errors.close()

final_variables = []
for df_heptads in os.listdir(os.path.join(os.getcwd(),"heptads_annotation")):
    df = pd.read_csv(f"heptads_annotation/{df_heptads}",sep="\t")
    variables = df.columns.tolist()
    for variable in variables:
        if variable not in final_variables:
            final_variables.append(variable)
final_variables = sorted(final_variables)
final_variables.remove("pdb")
final_variables.remove("Unnamed: 0")
final_variables.remove("net_charge")
final_variables.remove("unique_chains")

pdb_list = []
net_charge = []
unique_chains = []
df_output = pd.DataFrame()
for variable in final_variables:
    globals()[f"{variable}"] = []

for df_heptads in os.listdir(os.path.join(os.getcwd(),"heptads_annotation")):
    df = pd.read_csv(f"heptads_annotation/{df_heptads}", sep="\t")
    columns = df.columns.tolist()
    pdb_list.append(df.loc[0,"pdb"])
    net_charge.append(df.loc[0, "net_charge"])
    unique_chains.append(df.loc[0, "unique_chains"])
    for variable in final_variables:
        if variable in columns:
            globals()[f"{variable}"].append(df.loc[0,f"{variable}"])
        else:
            globals()[f"{variable}"].append(" ")
df_output["pdb"] = pdb_list
df_output["net_charge"] = net_charge
df_output["unique_chains"] = unique_chains
for variable in final_variables:
    df_output[f"{variable}"] = globals()[f"{variable}"]
print(df_output)

df_output.to_csv("output_descriptors.csv", sep="\t")