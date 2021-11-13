from parserData import *
import numpy as np

coordinates = []
elements = []

for line in pdb1s5l.split("\n"):
    if line.startswith("ATOM"):
        entries = line.split(" ")
        result = []
        for entry in entries:
            if entry != "":
                result.append(entry)
        (x, y, z, w) = (result[6], result[7], result[8], 1)
        coordinates.append((x,y,z, w))
        elements.append(result[-1])

coordinatesInNumpy, elementsInNumpy = np.array(coordinates), np.array(elements)
print(coordinatesInNumpy, elementsInNumpy)