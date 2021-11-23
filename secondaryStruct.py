# Name: Claire Chen
# Andrew ID: ccz
from cmu_112_graphics import *
from seqAnalysis import *
import numpy as np
import parser
from PIL import Image
import math

coordinatesInNumpy = parser.coordinatesInNumpy
elemsInNumpy = parser.elemsInNumpy
aminoAcidSeq = parser.aminoAcidSeq
helices = parser.helices
sheets = parser.sheets

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
    app.width, app.height = 1000, 800
    # Source of image: https://i.stack.imgur.com/T57R5.png (rescaled)
    app.image_helix = Image.open("helix.png")
    app.structHelix = []
    for helix in helices:
        startHelixPos, endHelixPos = helices[helix][1], helices[helix][3]
        if startHelixPos in aminoAcidSeq and endHelixPos in aminoAcidSeq:
            startAAPos = aminoAcidSeq[startHelixPos][0]
            endAAPos = aminoAcidSeq[endHelixPos][-1]
            startPos = threeDToTwoD(startAAPos)
            endPos = threeDToTwoD(endAAPos)
            app.structHelix.append((startPos, endPos))
    # Source of image: myself (using Adobe Sketch)
    app.image_sheet = Image.open("sheet.png")
    app.structSheet = []
    for sheet in sheets:
        startSheetPos, endSheetPos = sheets[sheet][1], sheets[sheet][3]
        if startSheetPos in aminoAcidSeq and endSheetPos in aminoAcidSeq:
            startAAPos = aminoAcidSeq[startSheetPos][0]
            endAAPos = aminoAcidSeq[endSheetPos][-1]
            startPos = threeDToTwoD(startAAPos)
            endPos = threeDToTwoD(endAAPos)
            app.structSheet.append((startPos, endPos))

def redrawAll(app, canvas):  
    for struct in app.structHelix:
        startPos, endPos = struct[0], struct[1]
        drawHelix(app, canvas, startPos, endPos)
    for struct in app.structSheet:
        startPos, endPos = struct[0], struct[1]
        drawSheet(app, canvas, startPos, endPos)
        
def drawHelix(app, canvas, startPos, endPos):
    x0, x1 = startPos[0] + app.width / 1.8, endPos[0] + app.width / 1.8
    y0, y1 = startPos[1] - app.height / 4, endPos[1] - app.height / 4
    # canvas.create_line(x0, y0, x1, y1, fill = "blue", width = 2)
    # The original image is 1600 pixels wide (we only care about the length)
    dx, dy = x1 - x0, y1 - y0
    distance = (dx ** 2 + dy ** 2) ** 0.5
    ratio = distance / 1600
    image = app.image_helix
    radOfRot = math.atan(dy/dx)
    degreeOfRot = math.degrees(radOfRot)
    imageOfHelix = image.rotate(degreeOfRot)
    # The ratio is further scaled since the original ratio would render 
    # the helices too large
    image_helix_scaled = app.scaleImage(imageOfHelix, ratio)
    # Anchor the center of the image to the midpoint between the start 
    # and end pos
    center = ((x0 + x1) / 2, (y0 + y1) / 2)
    canvas.create_image(center[0], center[1], 
                        image=ImageTk.PhotoImage(image_helix_scaled))

def drawSheet(app, canvas, startPos, endPos):
    x0, x1 = startPos[0] + app.width / 1.8, endPos[0] + app.width / 1.8
    y0, y1 = startPos[1] - app.height / 4, endPos[1] - app.height / 4
    # canvas.create_line(x0, y0, x1, y1, fill = "blue", width = 2)
    # The original image is 1998 pixels wide (we only care about the length)
    dx, dy = x1 - x0, y1 - y0
    distance = (dx ** 2 + dy ** 2) ** 0.5
    ratio = distance / 1998
    image = app.image_sheet
    radOfRot = math.atan(dy/dx)
    degreeOfRot = math.degrees(radOfRot)
    imageOfSheet = image.rotate(degreeOfRot)
    # The ratio is further scaled since the original ratio would render 
    # the helices too large
    image_sheet_scaled = app.scaleImage(imageOfSheet, ratio)
    # Anchor the center of the image to the midpoint between the start 
    # and end pos
    center = ((x0 + x1) / 2, (y0 + y1) / 2)
    canvas.create_image(center[0], center[1], 
                        image=ImageTk.PhotoImage(image_sheet_scaled))

runApp(width = 1000, height = 800)