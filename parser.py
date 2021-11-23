import parserData
import numpy as np
from Bio.PDB import *

# protein = '1s5l'
# structure = MMCIFParser().get_structure('1s5l', 's5/1s5l.cif')
# pdb1s5l = str(structure)
pdb1s5l = parserData.pdb1s5l

coordinates = []
elements = []
helices = dict()
#         {Index: ["{FirstAA}", startPos,  "{LastAA}",  endPos, len]}
#  e.g.   {1:     ["ASN",          68,         "VAL",    78,    10 ]}
sheets = dict()
aminoAcidSeq = dict()

for line in pdb1s5l.split("\n"):
    if line.startswith("TER"):
        break
    elif line.startswith("ATOM"):
        entries = line.split(" ")
        result = []
        for entry in entries:
            if entry != "":
                result.append(entry)
        # Adds to dictionaries of coordinates and elements
        [x, y, z, w] = [float(result[6]), float(result[7]), float(result[8]), 1]
        coordinates.append([x,y,z, w])
        elements.append(result[-1])
        # Adds to indexed dictionary of amino acids for modeling secondary 
        # structures
        aaIndex = result[5]
        if aaIndex in aminoAcidSeq:
            [x, y, z, w] = [float(result[6]), float(result[7]), float(result[8]), 1]
            aminoAcidSeq[aaIndex].append([x, y, z, w])
        else:
            aminoAcidSeq[aaIndex] = []
    elif line.startswith("HELIX"):
        entries = []
        for entry in line.split(" "):
            if entry != "" and entry != "HELIX":
                entries.append(entry)
        entries = entries[1:3] + entries[4:6] + [entries[7]] + [entries[-1]]
        helices[entries[0]] = entries[1:]
    elif line.startswith("SHEET"):
        pass
        # entries = line.split(" ")
        # print(entries)
        # for entry in entries:
        #     if entry != "" and len(entry) != 1 and entry != entries[1]:
        #         pass

# print(aminoAcidSeq)
print(helices)
# print(sheets)

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
    row[0] = row[0] - middleZ
    row[1] = row[1] - middleZ
coordinatesInNumpy, elemsInNumpy = np.array(coordinates), np.array(elements)