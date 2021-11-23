from cmu_112_graphics import *
import math
import numpy as np
import parser

# coordinates = [[0, 0, 0, 1], [0, 100, 0, 1], [0, 0, 100, 1],
#               [100, 0, 0, 100], [100, 100, 0, 1], [100, 0, 100, 1],
#               [0, 100, 100, 1], [100, 100, 100, 1]]
# coordinatesInNumpy = np.array(coordinates)
# elems = ['N', 'C', 'O', 'H', 'C', 'P', 'H', 'S']
# elemsInNumpy = np.array(elems)

coordinatesInNumpy = parser.coordinatesInNumpy
elemsInNumpy = parser.elemsInNumpy

class Protein(object):
    def __init__(self, coordinate, atom):
        location = [coordinate[0], coordinate[1], 
                    coordinate[2], coordinate[3]]
        self.coordinate = np.array(location)
        self.staticCoords = np.array(location)
        self.coordinate2D = threeDToTwoD(self.coordinate)
        self.atom = atom

    def scaleUp(self):
        entries = [[1.05, 0, 0, 0], 
                   [0, 1.05, 0, 0],
                   [0, 0, 1.05, 0],
                   [0, 0, 0, 1]]
        scaleUpMatrix = np.array(entries)
        product = np.matmul(scaleUpMatrix, self.coordinate)
        self.coordinate = product
    
    def scaleDown(self):
        entries = [[0.95, 0, 0, 0], 
                   [0, 0.95, 0, 0],
                   [0, 0, 0.95, 0],
                   [0, 0, 0, 1]]
        scaleUpMatrix = np.array(entries)
        product = np.matmul(scaleUpMatrix, self.coordinate)
        self.coordinate = product

    def rotateAroundX(self, sign):
        entries = [[math.cos(sign * math.pi / 12), -math.sin(sign * math.pi / 12), 0, 0],
                   [math.sin(sign * math.pi / 12), math.cos(sign * math.pi / 12), 0, 0],
                   [0, 0, 1, 0],
                   [0, 0, 0, 1]]
        rotateMatrix = np.array(entries)
        product = np.matmul(rotateMatrix, self.coordinate)
        self.coordinate = product
    def rotateAroundY(self, sign):
        entries = [[1, 0, 0, 0],
                   [0, math.cos(sign * math.pi / 12), math.sin(sign * math.pi / 12), 0],
                   [0, -math.sin(sign * math.pi / 12), math.cos(sign * math.pi / 12), 0],
                   [0, 0, 0, 1]]
        rotateMatrix = np.array(entries)
        product = np.matmul(rotateMatrix, self.coordinate)
        self.coordinate = product
    def rotateAroundZ(self, sign):
        entries = [[math.cos(sign * math.pi / 12), 0, math.sin(sign * math.pi / 12), 0],
                   [0, 1, 0, 0],
                   [-math.sin(sign * math.pi / 12), 0, math.cos(sign * math.pi / 12), 0],
                   [0, 0, 0, 1]]
        rotateMatrix = np.array(entries)
        product = np.matmul(rotateMatrix, self.coordinate)
        self.coordinate = product

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

def redrawAll(app, canvas):
    drawAxes(app, canvas)
    for atom in app.atoms:
        # coordinate = (atom.coordinate2D[0] + app.width / 2, 
        #               atom.coordinate2D[1] + 4 * app.height / 10)
        coordinate = (atom.coordinate2D[0] + app.width / 2, 
                      atom.coordinate2D[1])
        if atom.atom == 'C':
            drawCarbon(app, canvas, coordinate)
        elif atom.atom == 'H':
            drawHydrogen(app, canvas, coordinate)
        elif atom.atom == 'O':
            drawOxygen(app, canvas, coordinate)
        elif atom.atom == 'N':
            drawNitrogen(app, canvas, coordinate)
        elif atom.atom == 'P':
            drawPhosphorus(app, canvas, coordinate)
        elif atom.atom == 'S':
            drawSulfur(app, canvas, coordinate)
    canvas.create_text(app.width / 6, 5 * app.height / 6, text = f"app level = {app.level}")
def drawAxes(app, canvas):
    # X axis
    canvas.create_line(app.width / 5, 4 * app.height / 5, app.width / 2, app.height / 2)
    # Y axis
    canvas.create_line(app.width / 2, app.height / 2, 4 * app.width / 5, 4 * app.height / 5)
    # Z axis
    canvas.create_line(app.width / 2, app.height / 5, app.width / 2, app.height / 2)

def drawCarbon(app, canvas, coordinate):
    r = 2 * (1 + 0.001 * app.level)
    (x, y) = (coordinate[0], coordinate[1])
    canvas.create_oval(x-r, y-r, x+r, y+r, fill = 'black')

def drawHydrogen(app, canvas, coordinate):
    r = 1 * (1 + 0.001 * app.level)
    (x, y) = (coordinate[0], coordinate[1])
    canvas.create_oval(x-r, y-r, x+r, y+r, fill = 'floral white', outline = 'floral white')

def drawOxygen(app, canvas, coordinate):
    r = 2 * (1 + 0.001 * app.level)
    (x, y) = (coordinate[0], coordinate[1])
    canvas.create_oval(x-r, y-r, x+r, y+r, fill = 'IndianRed3', outline = "IndianRed3")

def drawNitrogen(app, canvas, coordinate):
    r = 2 * (1 + 0.001 * app.level)
    (x, y) = (coordinate[0], coordinate[1])
    canvas.create_oval(x-r, y-r, x+r, y+r, fill = 'RoyalBlue3', outline = 'RoyalBlue3')

def drawPhosphorus(app, canvas, coordinate):
    r = 2 * (1 + 0.001 * app.level)
    (x, y) = (coordinate[0], coordinate[1])
    canvas.create_oval(x-r, y-r, x+r, y+r, fill = 'gold', outline = "gold")

def drawSulfur(app, canvas, coordinate):
    r = 2 * (1 + 0.001 * app.level)
    (x, y) = (coordinate[0], coordinate[1])
    canvas.create_oval(x-r, y-r, x+r, y+r, fill = 'LightPink1', outline = 'LightPink1')