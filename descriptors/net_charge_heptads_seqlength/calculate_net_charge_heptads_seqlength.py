# !usr/bin/env python

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


# ARGUMENTS:
fullPDB = sys.argv[1]
list = sys.argv[2]
df = pd.read_csv(list, sep="\t")
dir = os.getcwd()
if not "fasta_files" in os.listdir(dir):
    subprocess.run(f"mkdir fasta_files/", shell=True)
if not "modified_PDBs" in os.listdir(dir):
    subprocess.run(f"mkdir modified_PDBs/", shell=True)
resultsHeptads = open("output_heptads_annotation.txt", "w")


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
def assign_heptads(pdb):
    heptads = str()
    #hydrophobicResidues = str()
    if not f"{pdb}.fasta" in os.listdir(f"{dir}/fasta_files/"):
        subprocess.run(f"wget -O ./fasta_files/{pdb}.fasta https://www.rcsb.org/fasta/entry/{pdb}", shell=True)
    fasta_file = f"{dir}/fasta_files/{pdb}.fasta"
    fasta_file = open(fasta_file, "r")
    fasta = fasta_file.read().splitlines()
    for line in fasta:
        i = 0
        if ">" not in line:
            seq = line
            seq_length = int(len(seq))
            print(seq_length)
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
            resultsHeptads.write(f"\n{pdb}\n{seq}\n{heptads}\n")
        # for ele in heptads:
        #     if ele == "a" or ele == "d":
        #         pos = heptads.index(ele)
        #         hydrophobicResidues = hydrophobicResidues+seq[pos]
            return seq, positions, seq_length#, hydrophobicResidues
    fasta_file.close()


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
            df.at[i, "net_charge"] = calculate_net_charge(ent_file)
            try: # Assign heptad repeats
                structure = parser.get_structure(pdb,ent_file)
                print(f"\n{pdb}")
                sequence, positions, seq_length = assign_heptads(pdb)
                #df.at[i, "sequence"] = sequence
                df.at[i, "sequence_length"] = seq_length
                if len(positions) > 0:
                    df.at[i, "heptads"] = "yes"
                else:
                    df.at[i, "heptads"] = "no"
                df.at[i, "heptad_number"] = len(positions)
                df.at[i, "heptad_positions"] = str(positions)
                ########################3# sr.compute(structure, level="S")
                # print(round(structure.sasa,2))
                # if hydrophobicCore:
                #     try: # Calculate hydrophobic core's net-charge
                #         df.at[i, "hydrophobic_core(HC)"] = hydrophobicCore
                #         df.at[i, "HC_net_charge"] = calculate_net_charge(hydrophobicCore)
                #     except:
                #         print(f"Error in {pdb}: impossible to calculate the net charge.")
                try:
                    heptadRepeats = bfactorToHeptadAnnotation(structure, pdb, positions, sequence)
                    df.at[i, "heptad_repeats"] = ",".join(heptadRepeats)
                except:
                    print(f"Error in {pdb}: impossible to modify b-factors.")
            except:
                print(f"Error in {pdb}: impossible to assign heptad repeats.")
        except:
            print("Error in {pdb}: impossible to calculate global net charge.")
    except:
        print(f"Error: {pdb} PDB not found!")


df.to_csv("output_descriptors.csv", sep="\t")

resultsHeptads.close()