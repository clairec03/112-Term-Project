from parserData import *
import numpy as np
from Bio.PDB import *
import os
import subprocess

# protein = '1s5l'
# structure = MMCIFParser().get_structure('1s5l', 's5/1s5l.cif')
# print(structure)
# # Run the bash script that converts the ENT file into a pure text file
# print("start")
# subprocess.call("toPDB.sh")
# print("end")

coordinates = []
elements = []

for line in pdb1s5l.split("\n"):
    if line.startswith("TER"):
        break
    elif line.startswith("ATOM"):
        entries = line.split(" ")
        result = []
        for entry in entries:
            if entry != "":
                result.append(entry)
        [x, y, z, w] = [float(result[6]), float(result[7]), float(result[8]), 1]
        coordinates.append([x,y,z, w])
        elements.append(result[-1])

testArray = np.array(coordinates)
# Finds the z-coordinate of the lowest and highest atom
minZ = 1000000000
maxZ = -1
for row in coordinates:
    currZ = row[2]
    if currZ < minZ:
        minZ = currZ
    if currZ > maxZ:
        maxZ = currZ
# Subtract the middle z value from all z-coordinates to ensure the protein is 
# shifted downwards to ensure a proper view of the molecule
middleZ = (maxZ + minZ) / 2
for row in coordinates:
    row[2] = row[2] - middleZ
coordinatesInNumpy, elemsInNumpy = np.array(coordinates), np.array(elements)