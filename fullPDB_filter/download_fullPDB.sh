#!/bin/bash

# Usage: sh script.sh
# Aim: to download the full PDB. Files are saved with .ent format.

nohup wget -r ftp://ftp.wwpdb.org/pub/pdb/data/structures/all/pdb/ &> log

# After everything is downloaded, execute the following command:
gunzip *.gz
