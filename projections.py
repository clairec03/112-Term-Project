from cmu_112_graphics import *
import math
import numpy as np
import parser


class Protein(object):
    def __init__(self, coordinate, atom):
        location = [coordinate[0], coordinate[1], 
                    coordinate[2], coordinate[3]]
        self.coordinate = np.array(location)
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

class Helix(object):
    def __init__(self, startingCoords, endingCoords):
        self.start = startingCoords
        self.end = endingCoords
        self.start2D = threeDToTwoD(self.start)
        self.end2D = threeDToTwoD(self.end)

    def scaleUp(self):
        entries = [[1.05, 0, 0, 0], 
                   [0, 1.05, 0, 0],
                   [0, 0, 1.05, 0],
                   [0, 0, 0, 1]]
        scaleUpMatrix = np.array(entries)
        productStart = np.matmul(scaleUpMatrix, self.start)
        productEnd = np.matmul(scaleUpMatrix, self.end)
        self.start = productStart
        self.end = productEnd

    def scaleDown(self):
        entries = [[0.95, 0, 0, 0], 
                   [0, 0.95, 0, 0],
                   [0, 0, 0.95, 0],
                   [0, 0, 0, 1]]
        scaleUpMatrix = np.array(entries)
        productStart = np.matmul(scaleUpMatrix, self.start)
        productEnd = np.matmul(scaleUpMatrix, self.end)
        self.start = productStart
        self.end = productEnd

    def rotateAroundX(self, sign):
        entries = [[math.cos(sign * math.pi / 12), -math.sin(sign * math.pi / 12), 0, 0],
                   [math.sin(sign * math.pi / 12), math.cos(sign * math.pi / 12), 0, 0],
                   [0, 0, 1, 0],
                   [0, 0, 0, 1]]
        rotateMatrix = np.array(entries)
        productStart = np.matmul(rotateMatrix, self.start)
        productEnd = np.matmul(rotateMatrix, self.end)
        self.start = productStart
        self.end = productEnd

    def rotateAroundY(self, sign):
        entries = [[1, 0, 0, 0],
                   [0, math.cos(sign * math.pi / 12), math.sin(sign * math.pi / 12), 0],
                   [0, -math.sin(sign * math.pi / 12), math.cos(sign * math.pi / 12), 0],
                   [0, 0, 0, 1]]
        rotateMatrix = np.array(entries)
        productStart = np.matmul(rotateMatrix, self.start)
        productEnd = np.matmul(rotateMatrix, self.end)
        self.start = productStart
        self.end = productEnd

    def rotateAroundZ(self, sign):
        entries = [[math.cos(sign * math.pi / 12), 0, math.sin(sign * math.pi / 12), 0],
                   [0, 1, 0, 0],
                   [-math.sin(sign * math.pi / 12), 0, math.cos(sign * math.pi / 12), 0],
                   [0, 0, 0, 1]]
        rotateMatrix = np.array(entries)
        productStart = np.matmul(rotateMatrix, self.start)
        productEnd = np.matmul(rotateMatrix, self.end)
        self.start = productStart
        self.end = productEnd

class Sheet(Helix):
    def __init__(self, startingCoords, endingCoords):
        super().__init__(startingCoords, endingCoords)
        self.start = startingCoords
        self.end = endingCoords
        self.start2D = threeDToTwoD(self.start)
        self.end2D = threeDToTwoD(self.end)

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