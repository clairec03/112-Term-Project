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

def appStarted(app):
    app.inputs = []
    app.atoms = list()
    app.coordinates, app.elems = coordinatesInNumpy, elemsInNumpy
    for i in range(len(coordinatesInNumpy)):
        app.atoms.append(Protein(app.coordinates[i], app.elems[i]))
    app.width, app.height = 1000, 800

def redrawAll(app, canvas):  
    for index in range(len(app.atoms)-1):
        coords1 = (5 * app.atoms[index].coordinate2D[0] + app.width / 1.5, 
                    5 * app.atoms[index].coordinate2D[1] + app.height / 2)
        coords2 = (5 * app.atoms[index+1].coordinate2D[0] + app.width / 1.5, 
                    5 * app.atoms[index+1].coordinate2D[1] + app.height / 2)
        drawHelix(canvas, coords1, coords2)

def drawHelix(canvas, coords1, coords2):
    print((coords1[0]))
    canvas.create_line(coords1[0], coords2[1], coords2[0], coords2[1],
                        width=1, fill = "blue")

runApp(width = 1000, height = 800)