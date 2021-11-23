# Name: Claire Chen
# Andrew ID: ccz
from cmu_112_graphics import *
from seqAnalysis import *
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
    # Then, do an orthographic projection on the product (with scaling)
    projections = [[8, 0, 0],
                  [0, 8, 0],
                  [0, 0, 0]]
    projectToXY = np.array(projections)
    result = np.matmul(projectToXY, product)
    return result[0:2]

def appStarted(app):
    app.inputs = []
    app.atoms = list()
    app.coordinates, app.elems = coordinatesInNumpy, elemsInNumpy
    app.width, app.height = 1000, 800
    app.struct = []
    # print(helices)
    for helix in helices:
        print(helices[helix][1])
        startHelixPos, endHelixPos = helices[helix][1], helices[helix][3]
        if startHelixPos in aminoAcidSeq and endHelixPos in aminoAcidSeq:
            startAAPos = aminoAcidSeq[startHelixPos][0]
            endAAPos = aminoAcidSeq[endHelixPos][-1]
            startPos = threeDToTwoD(startAAPos)
            endPos = threeDToTwoD(endAAPos)
            app.struct.append((startPos, endPos))
        # struct = []
        # for index in range(int(helices[helix][1]), int(helices[helix][3]) + 1):
        #     # print(index)
        #     if index in aminoAcidSeq:
        #         coords = aminoAcidSeq[index]
        #         aaCoords = threeDToTwoD(coords)
        #         # print(aaCoords)
        #         struct.append(aaCoords)
    url_helix = "https://upload.wikimedia.org/wikipedia/commons/thumb/6/65/A-Helix.svg/744px-A-Helix.svg.png"
    app.image_helix = app.loadImage(url_helix)

def redrawAll(app, canvas):  
    for struct in app.struct:
        startPos, endPos = struct[0], struct[1]
        drawHelix(app, canvas, startPos, endPos)
        
def drawHelix(app, canvas, startPos, endPos):
    x0, x1 = startPos[0] + app.width / 1.8, endPos[0] + app.width / 1.8
    y0, y1 = startPos[1] - app.height / 4, endPos[1] - app.height / 4
    # canvas.create_line(x0, y0, x1, y1, fill = "blue", width = 2)
    # The original image is 744 * 177 pixels (we only care about the length)
    distance = ((x0 - x1) ** 2 + (y0 - y1) ** 2) ** 0.5
    ratio = distance / 744
    image_helix_scaled = app.scaleImage(app.image_helix, ratio)
    # Anchor the center of the image to the midpoint between the start 
    # and end pos
    center = ((x0 + x1) / 2, (y0 + y1) / 2)
    canvas.create_image(center[0], center[1], 
                        image=ImageTk.PhotoImage(image_helix_scaled))

runApp(width = 1000, height = 800)