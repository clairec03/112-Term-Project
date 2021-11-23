# Name: Claire Chen
# Andrew ID: ccz
from cmu_112_graphics import *
from seqAnalysis import *
import numpy as np
from projections import *
import numpy as np
import parser

coordinatesInNumpy = parser.coordinatesInNumpy
elemsInNumpy = parser.elemsInNumpy
aminoAcidSeq = parser.aminoAcidSeq
helices = parser.helices

def threeDToTwoD(coordinate):
    # Matrix from https://en.wikipedia.org/wiki/Isometric_projection
    scalar = 1 / (6 ** 0.5)
    entries = [[scalar * 3 ** 0.5, 0, -scalar * 3 ** 0.5],
                [scalar * 1, scalar * 2, scalar *1],
                [scalar * 2 ** 0.5, -scalar * 2 ** 0.5, scalar * 2 ** 0.5]]
    matrix = np.array(entries)
    coordinates3D = coordinate[0:3]
    product = np.matmul(matrix, coordinates3D)
    # Then, do an orthographic projection on the product 
    projections = [[1, 0, 0],
                  [0, 1, 0],
                  [0, 0, 0]]
    projectToXY = np.array(projections)
    result = np.matmul(projectToXY, product)
    return result[0:2]

def appStarted(app):
    app.inputs = []
    app.atoms = list()
    app.coordinates, app.elems = coordinatesInNumpy, elemsInNumpy
    for i in range(len(coordinatesInNumpy)):
        app.atoms.append(Protein(app.coordinates[i], app.elems[i]))
    app.width, app.height = 1000, 800
    app.struct = []
    for helix in helices:
        struct = []
        for index in range(int(helices[helix][1]), int(helices[helix][3]) + 1):
            if index in aminoAcidSeq:
                coords = aminoAcidSeq[index]
                for coords in aminoAcidSeq[index]:
                    aaCoords = threeDToTwoD(coords)
                    struct.append(aaCoords)
        print(struct)

# def redrawAll(app, canvas):  
#         drawHelix(canvas, struct)
#         print(struct)
        
def drawHelix(canvas, struct):
    for point in struct:
        # print(point)
        x, y, r = point[0], point[1], 1
        canvas.create_oval(x - r, y - r, x + r, y + r, fill = "blue", 
                                                    outline = "blue")

runApp(width = 1000, height = 800)