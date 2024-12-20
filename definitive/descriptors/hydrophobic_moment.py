# !usr/bin/env python
import traceback

# Usage: python script.py fullPDB_path PDB_list(.txt)

# LIBRARIES IMPORTATION:
from Bio.PDB import *
parser = PDBParser(QUIET=True)
io = PDBIO()
sr = ShrakeRupley()
import os, sys, subprocess
import pandas as pd
import re
import sys
import matplotlib.pyplot as plt
import random
from scipy.fft import fft, fftfreq
import numpy as np

# ARGUMENTS:
fullPDB = sys.argv[1]
df = sys.argv[2]
df = pd.read_csv(df, sep="\t")
dir = os.getcwd()


# DICTIONARIES (HYDROPHOBICITY SCALES):
# We create dictionaries to store the hydropathy indexes of the residues based on different hydrophobicity scales:
kyte_doolittle = {
	"ALA": 1.8,
	"CYS": 2.5,
	"ASP": -3.5,
	"GLU": -3.5,
	"PHE": 2.8,
	"GLY": -0.4,
	"HIS": -3.2,
	"ILE": 4.5,
	"LYS": -3.9,
	"LEU": 3.8,
	"MET": 1.9,
	"ASN": -3.5,
	"PRO": -1.6,
	"GLN": -3.5,
	"ARG": -4.5,
	"SER": -0.8,
	"THR": -0.7,
	"VAL": 4.2,
	"TRP": -0.9,
	"TYR": -1.3
}

hopp_woods = {
	"ALA": -0.5,
	"CYS": -1.0,
	"ASP": 3.0,
	"GLU": 3.0,
	"PHE": -2.5,
	"GLY": 0.0,
	"HIS": -0.5,
	"ILE": -1.8,
	"LYS": 3.0,
	"LEU": -1.8,
	"MET": -1.3,
	"ASN": 0.2,
	"PRO": 0.0,
	"GLN": 0.2,
	"ARG": 3.0,
	"SER": 0.3,
	"THR": -0.4,
	"VAL": -1.5,
	"TRP": -3.4,
	"TYR": -2.3,
}

cornette = {
	"ALA": 0.2,
	"CYS": 4.1,
	"ASP": -3.1,
	"GLU": -1.8,
	"PHE": 4.4,
	"GLY": 0.0,
	"HIS": 0.5,
	"ILE": 4.8,
	"LYS": -3.1,
	"LEU": 5.7,
	"MET": 4.2,
	"ASN": -0.5,
	"PRO": -2.2,
	"GLN": -2.8,
	"ARG": 1.4,
	"SER": -0.5,
	"THR": -1.9,
	"VAL": 4.7,
	"TRP": 1.0,
	"TYR": 3.2
}

eisenberg = {
	"ALA": 0.62,
	"CYS": 0.29,
	"ASP": -0.90,
	"GLU": -0.74,
	"PHE": 1.19,
	"GLY": 0.48,
	"HIS": -0.40,
	"ILE": 1.38,
	"LYS": -1.50,
	"LEU": 1.06,
	"MET": 0.64,
	"ASN": -0.78,
	"PRO": 0.12,
	"GLN": -0.85,
	"ARG": -2.53,
	"SER": -0.18,
	"THR": -0.05,
	"VAL": 1.08,
	"TRP": 0.81,
	"TYR": 0.26
}

rose = {
	"ALA": 0.74,
	"CYS": 0.91,
	"ASP": 0.62,
	"GLU": 0.62,
	"PHE": 0.88,
	"GLY": 0.72,
	"HIS": 0.78,
	"ILE": 0.88,
	"LYS": 0.52,
	"LEU": 0.85,
	"MET": 0.85,
	"ASN": 0.63,
	"PRO": 0.64,
	"GLN": 0.62,
	"ARG": 0.64,
	"SER": 0.66,
	"THR": 0.70,
	"VAL": 0.86,
	"TRP": 0.85,
	"TYR": 0.76
}

janin = {
	"ALA": 0.3,
	"CYS": 0.9,
	"ASP": -0.6,
	"GLU": -0.7,
	"PHE": 0.5,
	"GLY": 0.3,
	"HIS": -0.1,
	"ILE": 0.7,
	"LYS": -1.8,
	"LEU": 0.5,
	"MET": 0.4,
	"ASN": -0.5,
	"PRO": -0.3,
	"GLN": -0.7,
	"ARG": -1.4,
	"SER": -0.1,
	"THR": -0.2,
	"VAL": 0.6,
	"TRP": 0.3,
	"TYR": -0.4
}

engelman_ges = {
	"ALA": 1.6,
	"CYS": 2.0,
	"ASP": -9.2,
	"GLU": -8.2,
	"PHE": 3.7,
	"GLY": 1.0,
	"HIS": -3.0,
	"ILE": 3.1,
	"LYS": -8.8,
	"LEU": 2.8,
	"MET": 3.4,
	"ASN": -4.8,
	"PRO": -0.2,
	"GLN": -4.1,
	"ARG": -12.3,
	"SER": 0.6,
	"THR": 1.2,
	"VAL": 2.6,
	"TRP": 1.9,
	"TYR": -0.7
}

consensus1982 = {
	"ALA": 0.25,
	"CYS": 0.04,
	"ASP": -0.72,
	"GLU": -0.62,
	"PHE": 0.61,
	"GLY": 0.16,
	"HIS": -0.40,
	"ILE": 0.73,
	"LYS": -1.10,
	"LEU": 0.53,
	"MET": 0.26,
	"ASN": -0.64,
	"PRO": -0.07,
	"GLN": -0.69,
	"ARG": -1.76, # in the original paper the value is rounded to -1.8
	"SER": -0.26,
	"THR": -0.18,
	"VAL": 0.54,
	"TRP": 0.37,
	"TYR": 0.02
}

aminoacids = [
    "ALA",
    "ARG",
    "ASN",
    "ASP",
    "CYS",
    "GLN",
    "GLU",
    "GLY",
    "HIS",
    "ILE",
    "LEU",
    "LYS",
    "MET",
    "PHE",
    "PRO",
    "SER",
    "THR",
    "TRP",
    "TYR",
    "VAL",
    "MSE"
]

aminoacids_short = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V", "MSE": "U"
}

# Hydrophobic (Nonpolar) residues
hydrophobic = ["A", "V", "L", "I", "F", "W", "M", "P","C","G","U"]
# Charged residues (Basic and Acidic)
charged = ["K", "R", "H","D", "E"]
# Aromatic residues
aromatic = ["F", "W", "Y"]
# Polar residues
polar = ["S", "T", "Q", "N"]


# Diccionario one-hot para aminoácidos
aminoacids1 = ["A", "R", "N", "D", "C", "Q", "E", "G", "H", "I","L", "K", "M", "F", "P", "S", "T", "W", "Y", "V", "U"]
one_hot_dict = {aa: np.array([1 if i == idx else 0 for i in range(len(aminoacids1))]) for idx, aa in enumerate(aminoacids1)}

# FUNCTIONS:
# Función para codificar una secuencia de aminoácidos
def encode_sequence(seq):
	encoded = np.array([one_hot_dict[aa] for aa in seq])
	# Calcular el promedio de los vectores one-hot
	averaged_vector = np.mean(encoded, axis=0)
	return averaged_vector



'''
ADD GRAPHICAL REPRESENTATION AS SEEN IN PAPER
'''
# def calculate_hydrophobic_moment(ent):
# 	# Amphiphilic (amphiphatic): A molecule having both hydrophobic (nonpolar) and hydrophilic (polar) regions.
# 	return "aaaaa"
# def fourier_transform():
# 	# Number of sample points
# 	N = 20
# 	# sample spacing
# 	T = 4 / 10
# 	x = np.linspace(0.0, N * T, N, endpoint=False)
# 	y = np.sin(50.0 * 2.0 * np.pi * x) + 0.5 * np.sin(80.0 * 2.0 * np.pi * x)
# 	yf = fft(y)
# 	xf = fftfreq(N, T)[:N // 2]
# 	plt.plot(xf, 2.0 / N * np.abs(yf[0:N // 2]))
# 	plt.grid()
# 	plt.show()

# EXECTION:
# Iteration through all the PDBs list:
counter = 0
errors_num = 0
for i in df.index:
	pdb = df.loc[i,"pdb"]
	counter += 1
	residues_list = []
	try:
		ent_file = os.path.join(fullPDB,f"pdb{pdb}.ent")
		structure = parser.get_structure(f"{pdb}", f"{ent_file}")
		for model in structure:
			for chain in model:
				if Chain.Chain.get_id(chain) == "A":
					for residue in chain:
						res_name = residue.get_resname()
						if res_name in aminoacids:
							residues_list.append(res_name)
		# Constants
		angle_per_residue = 100 * (np.pi / 180)  # Convert to radians; 100 because we are dealing with helices; 160-180 for beta-sheet
		# Initialize x and y components
		x_sum = 0
		y_sum = 0
		# Loop through sequence
		"""
		AÑADIR FORMULA COGIENDO COORDENADAS DEL ARTÍCULO DE HYDROPHOBIC MOMENT,:
		https://github.com/JoaoRodrigues/hydrophobic_moment/blob/main/hydrophobic_moment.py -- github referencia
		CONSIDERAR DIRECCION VECTOR, DIFERENCIAS SEGÚN AREA INTERACCION CON MAS O MENOS HELICES
		
		- HELICAL PERFECTION: curvatura, más "ideales", geometria
		
		"""
		for index_res,residue in enumerate(residues_list):
			if residue in consensus1982:
				H_i = consensus1982[residue]
				theta_i = index_res * angle_per_residue
				x_sum += H_i * np.cos(theta_i)
				y_sum += H_i * np.sin(theta_i)
		# Compute hydrophobic moment
		hydrophobic_moment = round(np.sqrt(x_sum ** 2 + y_sum ** 2),2)
		print(hydrophobic_moment)
		df.at[i,"hydrophobic_moment"] = hydrophobic_moment

		# SEQUENCE AND HEPTADS ANNOTATION:
		seq = str()
		hydrophobic_core = 0
		prolines = 0
		heptads_core = 0
		aromatic_core = 0
		charged_core = 0
		for resi in residues_list:
			seq = seq + aminoacids_short[resi]
			if aminoacids_short[resi] in hydrophobic:
				hydrophobic_core += 1
			if aminoacids_short[resi] == "P":
				prolines += 1
			if aminoacids_short[resi] in aromatic:
				aromatic_core += 1
			if aminoacids_short[resi] in charged:
				charged_core += 1

		limit = 0
		heptads = 0
		max_heptads = round(len(seq)/7)
		#print("hi", max_heptads)
		#print(max_heptads)
		while limit < len(seq)-7:
			initial = seq[limit]
			last = seq[limit+3]
			if initial in hydrophobic and last in hydrophobic:
				heptads += 1
				limit += 7
			else:
				limit += 1
		# Embedding of sequence:
		df.at[i,"seq_length"] = len(residues_list)
		df.at[i, "seq"] = seq
		print(seq)
		# Lista de aminoácidos
		seq_encoded = encode_sequence(seq)
		print(seq_encoded)
		for index,ele in enumerate(seq_encoded):
			if ele != 0:
				df.at[i, f"res{index}"] = round(ele,2)
			else:
				df.at[i,f"res{index}"] = 0
		df.at[i,f"hydrophobic_res"] = round(hydrophobic_core/len(seq),2)
		df.at[i, f"disruptive_prolines"] = round(prolines/len(seq), 2)
		df.at[i, f"aromatic_res"] = round(aromatic_core / len(seq), 2)
		df.at[i, f"charged_res"] = round(charged_core / len(seq), 2)
		df.at[i, "canonical_heptads"] = round(heptads/max_heptads,2)
		print(heptads)
		print(f"{pdb.upper()} was correctly annotated. Stucture number {counter}.")
	except:
		errors_num += 1
		print(f"Error in {pdb.upper()}. Stucture number {counter}.")
df.to_csv(f"output_descriptors.csv",sep="\t")