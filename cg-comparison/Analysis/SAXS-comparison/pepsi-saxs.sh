#!/bin/sh

# Download PEPSI-SAXS from https://team.inria.fr/nano-d/software/pepsi-saxs/ and install on a system.
# It can be used for SAXS predictions as given below:

#./Pepsi-SAXS  <input PDB(s)> <experimental curve> [-o <output file>]
echo ./Pepsi-SAXS *.pdb SASDAC2.dat 
