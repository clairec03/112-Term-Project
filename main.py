# Name: Claire Chen
# Andrew ID: ccz
from cmu_112_graphics import *
import numpy as np
import math
from functools import cmp_to_key
from PIL import Image
from Bio.PDB import *
from pathlib import Path
import random

def appStarted(app):
    app.pdb = ''
    app.isParsing = False
    app.finishedParsing = False
    app.inputs = []
    app.input_dna = ''
    app.atoms = list()
    app.level = 50
    app.viewerCoords = (500, 400, 500)
    app.coordinates, app.elems = [], []
    app.aminoAcidSeq, app.helices, app.sheets = [], [], []
    # Source: https://www.biospace.com/getasset/aa8b842c-373c-4418-a466-a203bbe456d7/
    app.image_home = Image.open("homepage.jpg")
    app.image_home_scaled = app.scaleImage(app.image_home, 2/3)
    # Source: https://www.flaticon.com/free-icon/scissors_2168806?related_id=2168905&origin=search
    image_design_scissors = Image.open("scissors.png")
    app.img_scissors = app.scaleImage(image_design_scissors, 1/10)
    # Source: https://www.flaticon.com/free-icon/scissors_2168905?related_id=2168905&origin=search
    image_scissors_on_click = Image.open("scissors_outline.png")
    app.img_scissors_outline = app.scaleImage(image_scissors_on_click, 1/10)
    # Source: http://clipart-library.com/free/pencil-clipart-transparent-background.html
    image_design_pencil = Image.open("pencil.png")
    app.img_pencil = app.scaleImage(image_design_pencil, 1/10)
    app.buttonCenter1 = (2 * app.width / 5, 4.15 * app.height / 5,
                         3 * app.width / 5,  4.4 * app.height / 5)
    app.buttonTopLeft = (app.width / 12, app.height / 12, 
                         app.width / 4, app.height / 8)
    app.buttonBottomLeft = (app.width / 5, 3 * app.height / 4, 2 * app.width / 5, 
                            3.2 * app.height/4)
    app.buttonBottomRight = (3 * app.width/5, 3 * app.height / 4, 4 * app.width/5, 
                            3.2 * app.height/4)
    # Working on this rn
    app.buttonBottomMiddle = (app.width / 5, 3.3 * app.height / 4, 
                              2 * app.width / 5, 3.5 * app.height / 4)
    app.buttonSettings = (app.width / 10, app.height / 8, 
                          1.1 * app.width / 10, 1.08 * app.height / 8,)
    app.buttonInput = (app.width / 6, 7 * app.height / 8,
                            1.2 * app.width / 6, 7.5 * app.height / 8)
    app.buttonScissors = (5.2 * app.width / 6 - 25.6, 5.2 * app.height / 6 - 25.6, 
                        5.2 * app.width / 6 + 25.6, 5.2 * app.height / 6 + 25.6)
    app.buttonPencil = (4.8 * app.width / 6 - 25.6, 5.2 * app.height / 6 - 25.6, 
                        4.8 * app.width / 6 + 25.6, 5.2 * app.height / 6 + 25.6)
    app.dna_bases = {"adenine", "guanine", "cytosine", "thymine"}
    app.rna_bases = {"adenine", "guanine", "cytosine", "uracil"}
    app.amino_acids = {"START": {"AUG"},
                       "ALA": {"GCU", "GCC", "GCA", "GCG"},
                       "ARG": {"CGU", "CGC", "CGA", "CGG"},
                       "ASN": {"AAU", "AAC"},
                       "ASP": {"GAU", "GAC"},
                       "CYS": {"UGU", "UGC"},
                       "GLN": {"CAA", "CAG"},
                       "GLU": {"GAA", "GAG"},
                       "GLY": {"GGU", "GGC", "GGA", "GGG"},
                       "HIS": {"CAU", "CAC"},
                       "ILE": {"AUU", "AUC"},
                       "LEU": {"UUA", "UUG"},
                       "LYS": {"AAA", "AAG"},
                       "MET": {"AUG"},
                       "PHE": {"UUU", "UUC"},
                       "PRO": {"CCU", "CCC", "CCA", "CCG"},
                       "SER": {"AGU", "AGC"},
                       "THR": {"ACU", "ACC", "ACA", "ACG"},
                       "TRP": {"UGG"},
                       "TYR": {"UAU", "UAC"},
                       "VAL": {"GUU", "GUC", "GUA", "GUG"},
                       "STOP": {"UAA", "UAG"}
                       }
    app.width, app.height = 1000, 800
    app.wantInput = False
    app.timesteps = np.arange(0, 1000, .3)
    app.buttonX = (6 * app.width / 7, 2 * app.height / 5, 
                            6.5 * app.width / 7, 2.2 * app.height / 5)
    app.buttonY = (6 * app.width / 7, 2.3 * app.height / 5, 
                            6.5 * app.width / 7, 2.5 * app.height / 5)
    app.buttonZ = (6 * app.width / 7, 2.6 * app.height / 5, 
                            6.5 * app.width / 7, 2.8 * app.height / 5)
    app.buttonSwitch = (6 * app.width / 7, 3 * app.height / 5, 
                            6.5 * app.width / 7, 3.2 * app.height / 5)
    app.buttonXColor = "DeepSkyBlue2"
    app.buttonYColor = "DeepSkyBlue2"
    app.buttonZColor = "DeepSkyBlue2"
    app.upButtonColor = "SteelBlue1"
    app.downButtonColor = "SteelBlue1"
    app.leftButtonColor = "SteelBlue1"
    app.rightButtonColor = "SteelBlue1"
    app.buttonSwitchColor = "RosyBrown2"
    app.sheetColors = ("indian red", "pink", "navajo white", "cornflower blue",
                       "LightBlue1", "bisque", "PaleGreen1", "light slate blue", "thistle")
    app.leftArrowCoords = (6.9 * app.width / 8 - 30, 7 * app.height / 8 - 10,
                           6.9 * app.width / 8, 7 * app.height / 8 + 30)
    app.rightArrowCoords = (7.1 * app.width / 8 + 20, 7 * app.height / 8 - 10,
                            7.1 * app.width / 8 + 50, 7 * app.height / 8 + 30)
    app.downArrowCoords = (7 * app.width / 8 - 10, 7.3 * app.height / 8,
                           7 * app.width / 8 + 30, 7.3 * app.height / 8 + 30)
    app.upArrowCoords = (7 * app.width / 8 - 10, 6.9 * app.height / 8 - 30,
                         7 * app.width / 8 + 30, 6.9 * app.height / 8)
    app.isRotatingX, app.isRotatingY, app.isRotatingZ = False, False, False
    # Homepage buttons
    app.bottomLeftButtonColor = "steel blue"
    app.bottomRightButtonColor = "pale violet red"
    app.center1ButtonColor = "RoyalBlue4"
    app.center1ButtonMoved = False
    # Source of image: 
    # https://media.istockphoto.com/vectors/icon-flat-vector-id860692590?k=20&m=860692590&s=612x612&w=0&h=3xUcLjLCbCxrwMT6WGWhJCEVhDS6Znz0yNMslyywN0A=
    app.design = Image.open("helix_cover.jpg")
    app.image_design = app.scaleImage(app.design, 2/9)
    # Source of image:
    # https://www.science.org/do/10.1126/science.aal0634/abs/Fig_3_16x9.jpg
    app.protein = Image.open("protein_Model.jpg")
    app.protein_model = app.scaleImage(app.protein, 1/5)
    # Source of image: https://i.stack.imgur.com/T57R5.png (rescaled)
    app.image_helix = Image.open("helix.png")
    app.structHelix = []
    # Source of image: myself (using Adobe Sketch)
    app.image_sheet = Image.open("sheet.png")
    app.structSheet = []
    app.isAtomicModel = True
    app.isSecondaryStruct = False
    app.timerDelay = 1
    app.dna, app.rna, app.aminoAcids = '', '', ''
    app.blockZoomScale = 20
    app.spacing = app.blockZoomScale / 4
    app.y0 = 3 * app.height / 5
    app.y2 = 3 * app.height / 4
    app.y1 = (app.y0 + app.y2) / 2
    app.base = None
    # Source: https://www.cs.cmu.edu/~112/notes/notes-animations-part4.html#sidescrollerExamples
    app.scrollX = 0
    app.isDrawingHighlighterBox = False
    app.currIndex = 0
    app.missense_mutation, app.nonsense_mutation, app.frameshift = False, False, False
    resetApp(app)

def timerFired(app):
    if (app.isVisualPage and app.isAtomicModel and not app.isRotatingX and 
        not app.isRotatingY and not app.isRotatingZ):
        for atom in app.atoms:
            atom.rotateAroundZ(-.3)
            atom.coordinate2D = threeDToTwoD(atom.coordinate)
        # Source: https://stackoverflow.com/questions/5213033/sort-a-list-of-lists-with-a-custom-compare-function
        app.atoms = sorted(app.atoms, key=cmp_to_key(compare))
    elif (app.isVisualPage and app.isSecondaryStruct and not app.isRotatingX and 
        not app.isRotatingY and not app.isRotatingZ):
        for helix in app.structHelix:
            helix.rotateAroundZ(-1)
        for struct in app.structSheet:
            struct.rotateAroundZ(-1)
            
# Custom comparator function
def compare(atom1, atom2):
    return (getNorm(atom1.coordinate[0], atom1.coordinate[1], atom1.coordinate[2]) 
        - getNorm(atom2.coordinate[0], atom2.coordinate[1], atom2.coordinate[2]))

def mouseMoved(app, event):
    if app.buttonCenter1[0] <= event.x <= app.buttonCenter1[2]:
        app.center1ButtonColor = "midnight blue"
        app.center1ButtonMoved = True
    if inButton(event, app.buttonCenter1):
        app.center1ButtonColor = "midnight blue"
    if inButton(event, app.buttonTopLeft):
        pass
    if inButton(event, app.buttonTopLeft):
        pass
    if inButton(event, app.buttonBottomLeft):
        app.bottomLeftButtonColor = "midnight blue"
    if inButton(event, app.buttonBottomRight):
        app.bottomRightButtonColor = "IndianRed4"
    if inButton(event, app.buttonTopLeft):
        pass
    if inButton(event, app.buttonX):
        pass
    if inButton(event, app.buttonY):
        pass
    if inButton(event, app.buttonZ):
        pass
    if inButton(event, app.leftArrowCoords):
        pass
    if inButton(event, app.rightArrowCoords):
        pass
    if inButton(event, app.upArrowCoords):
        pass
    if inButton(event, app.downArrowCoords):
        pass
    else:
        app.bottomLeftButtonColor = "steel blue"
        app.bottomRightButtonColor = "pale violet red"
        app.center1ButtonColor = "RoyalBlue4"

def resetApp(app):
    app.isHomepage = True
    app.isDesignPage = False
    app.isHelpPage = False
    app.isIntroPage = False 
    app.drawVisualPage = False
    app.isVisualPage = False
    app.scissors_selected = False
    app.pencil_selected = False

def keyPressed(app, event):
    print("key pressed!")
    if event.key == "+":
        for atom in app.atoms:
            atom.scaleUp()
            atom.coordinate2D = threeDToTwoD(atom.coordinate)
        app.level += 1
    elif event.key == "-":
        for atom in app.atoms:
            atom.scaleDown()
            atom.coordinate2D = threeDToTwoD(atom.coordinate)
        app.level -= 1
    elif app.isDesignPage:
        if event.key == "Left":
                app.scrollX += 100
                app.changedBlocks = True
        elif event.key == "Right":
                app.scrollX -= 100
                app.changedBlocks = True
        elif event.key == "+":
            app.blockZoomScale *= 1.1
        elif event.key == "-":
            app.blockZoomScale *= 0.9
    elif app.isVisualPage:
        if event.key == "Left":
            if app.isRotatingX:
                if app.isAtomicModel:
                    for atom in app.atoms:
                        atom.rotateAroundX(-1)
                        atom.coordinate2D = threeDToTwoD(atom.coordinate)
                    app.atoms = sorted(app.atoms, key=cmp_to_key(compare))
                elif app.isSecondaryStruct:
                    for struct in app.structHelix:
                        struct.rotateAroundX(-1)
                        struct.start2D = threeDToTwoD(struct.start)
                        struct.end2D = threeDToTwoD(struct.end)
                    for struct in app.structSheet:
                        struct.rotateAroundX(-1)
            elif app.isRotatingY:
                if app.isAtomicModel:
                    for atom in app.atoms:
                        atom.rotateAroundY(-1)
                        atom.coordinate2D = threeDToTwoD(atom.coordinate)
                    app.atoms = sorted(app.atoms, key=cmp_to_key(compare))
                elif app.isSecondaryStruct:
                    for struct in app.structHelix:
                        struct.rotateAroundY(-1)
                        struct.start2D = threeDToTwoD(struct.start)
                        struct.end2D = threeDToTwoD(struct.end)
                    for struct in app.structSheet:
                        struct.rotateAroundY(-1)
            elif app.isRotatingZ:
                if app.isAtomicModel:
                    for atom in app.atoms:
                        atom.rotateAroundZ(-1)
                        atom.coordinate2D = threeDToTwoD(atom.coordinate)
                    app.atoms = sorted(app.atoms, key=cmp_to_key(compare))
                elif app.isSecondaryStruct:
                    for struct in app.structHelix:
                        struct.rotateAroundZ(-1)
                        struct.start2D = threeDToTwoD(struct.start)
                        struct.end2D = threeDToTwoD(struct.end)
                    for struct in app.structSheet:
                        struct.rotateAroundZ(-1)
            else:
                # Need to add the secondary structure rotations here
                for atom in app.atoms:
                    atom.rotateAroundZ(-1)
                    atom.coordinate2D = threeDToTwoD(atom.coordinate)
                    app.upButtonColor = "SteelBlue1"
                    app.downButtonColor = "SteelBlue1"
                    app.leftButtonColor = "SteelBlue4"
                    app.rightButtonColor = "SteelBlue1"
                app.atoms = sorted(app.atoms, key=cmp_to_key(compare))
        elif event.key == "Right":
            if app.isRotatingX:
                if app.isAtomicModel:
                    for atom in app.atoms:
                        atom.rotateAroundX(1)
                        atom.coordinate2D = threeDToTwoD(atom.coordinate)
                    app.atoms = sorted(app.atoms, key=cmp_to_key(compare))
                elif app.isSecondaryStruct:
                    for struct in app.structHelix:
                        struct.rotateAroundX(1)
                        struct.start2D = threeDToTwoD(struct.start)
                        struct.end2D = threeDToTwoD(struct.end)
                    for struct in app.structSheet:
                        struct.rotateAroundX(1)
            elif app.isRotatingY:
                if app.isAtomicModel:
                    for atom in app.atoms:
                        atom.rotateAroundY(1)
                        atom.coordinate2D = threeDToTwoD(atom.coordinate)
                    app.atoms = sorted(app.atoms, key=cmp_to_key(compare))
                elif app.isSecondaryStruct:
                    for struct in app.structHelix:
                        struct.rotateAroundY(1)
                        struct.start2D = threeDToTwoD(struct.start)
                        struct.end2D = threeDToTwoD(struct.end)
                    for struct in app.structSheet:
                        struct.rotateAroundY(1)
            elif app.isRotatingZ:
                if app.isAtomicModel:
                    for atom in app.atoms:
                        atom.rotateAroundZ(1)
                        atom.coordinate2D = threeDToTwoD(atom.coordinate)
                    app.atoms = sorted(app.atoms, key=cmp_to_key(compare))
                elif app.isSecondaryStruct:
                    for struct in app.structHelix:
                        struct.rotateAroundZ(1)
                        struct.start2D = threeDToTwoD(struct.start)
                        struct.end2D = threeDToTwoD(struct.end)
                    for struct in app.structSheet:
                        struct.rotateAroundZ(1)
            else:
                if app.isAtomicModel:
                    for atom in app.atoms:
                        atom.rotateAroundZ(1)
                        atom.coordinate2D = threeDToTwoD(atom.coordinate)
                        app.upButtonColor = "SteelBlue1"
                        app.downButtonColor = "SteelBlue1"
                        app.leftButtonColor = "SteelBlue1"
                        app.rightButtonColor = "SteelBlue4"
                    app.atoms = sorted(app.atoms, key=cmp_to_key(compare))
                elif app.isSecondaryStruct:
                    for struct in app.structHelix:
                        struct.rotateAroundZ(1)
                        struct.start2D = threeDToTwoD(struct.start)
                        struct.end2D = threeDToTwoD(struct.end)
                    for struct in app.structSheet:
                        struct.rotateAroundZ(1)
        elif event.key == "h":
            app.isHelpPage = True
        elif event.key == "Up":
            if app.isAtomicModel:
                for atom in app.atoms:
                    atom.rotateAroundX(-1)
                    atom.rotateAroundY(-1)
                    atom.coordinate2D = threeDToTwoD(atom.coordinate)
                    app.upButtonColor = "SteelBlue4"
                    app.downButtonColor = "SteelBlue1"
                    app.leftButtonColor = "SteelBlue1"
                    app.rightButtonColor = "SteelBlue1"
                app.atoms = sorted(app.atoms, key=cmp_to_key(compare))
            elif app.isSecondaryStruct:
                for struct in app.structHelix:
                    struct.rotateAroundX(-1)
                    struct.rotateAroundY(-1)
                    struct.start2D = threeDToTwoD(struct.start)
                    struct.end2D = threeDToTwoD(struct.end)
                for struct in app.structSheet:
                    struct.rotateAroundX(-1)
                    struct.rotateAroundY(-1)
        elif event.key == "Down":
            if app.isAtomicModel:
                for atom in app.atoms:
                    atom.rotateAroundX(1)
                    atom.rotateAroundY(1)
                    atom.coordinate2D = threeDToTwoD(atom.coordinate)
                    app.upButtonColor = "SteelBlue1"
                    app.downButtonColor = "SteelBlue4"
                    app.leftButtonColor = "SteelBlue1"
                    app.rightButtonColor = "SteelBlue1"
                app.atoms = sorted(app.atoms, key=cmp_to_key(compare))
            elif app.isSecondaryStruct:
                for struct in app.structHelix:
                    struct.rotateAroundX(1)
                    struct.rotateAroundY(1)
                    struct.start2D = threeDToTwoD(struct.start)
                    struct.end2D = threeDToTwoD(struct.end)
                for struct in app.structSheet:
                    struct.rotateAroundX(1)
    elif event.key == "h":
        app.isHelpPage = True
    elif app.isHelpPage and event.key == "x":
        app.isHelpPage = False
    elif app.isHomepage:
        app.isHomepage = False
        app.isIntroPage = True

def getUserInput(app, prompt):
    if app.wantInput:
        if app.isIntroPage:
            app.inputs.append(app.getUserInput(prompt).lower())
            if len(app.inputs) > 0 and len(app.inputs[-1]) and app.inputs[-1].isalnum():
                app.pdb = app.inputs[-1]
                processUserInput(app)
                app.isParsing = True
                app.wantInput = False
    if app.pencil_selected:
        input = app.getUserInput(prompt).upper()
        userInputToDNA(app, input)
        app.isDrawingHighlighterBox = False

def mousePressed(app, event):
    if inButton(event, app.buttonTopLeft):
        resetApp(app)
    if app.isHomepage:
        if inButton(event, app.buttonCenter1):
            app.isHomepage = False
            app.isDesignPage = True
    elif app.isDesignPage:
        if inButton(event, app.buttonTopLeft):
            resetApp(app)
        elif inButton(event, app.buttonScissors):
            app.scissors_selected = not app.scissors_selected
        elif inButton(event, app.buttonPencil):
            app.pencil_selected = not app.pencil_selected
        elif app.scissors_selected and event.y >= app.y0 and event.y <= app.y2:
            index = getIndex(app, event.x)
            app.dna = app.dna[0:index] + app.dna[index+1:-1]
            app.rna = transcribe(app.dna)
            app.aminoAcids = translate(app, app.rna)
            app.scissors_selected = False
        elif app.pencil_selected and event.y >= app.y0 and event.y <= app.y2:
            index = getIndex(app, event.x)
            app.currIndex = index
            app.isDrawingHighlighterBox = True
            getUserInput(app, "Please enter the DNA sequence to be inserted here\nWrite a string with A's, T's, G's, and C's only")
            app.dna = app.dna[0:index] + app.input_dna + app.dna[index:-1]
            getDNAComplement(app)
            app.rna = transcribe(app.dna)
            app.aminoAcids = translate(app, app.rna)
            app.pencil_selected = False
    elif app.isIntroPage:
        if inButton(event, app.buttonTopLeft):
            resetApp(app)
        elif inButton(event, app.buttonBottomLeft):
            if app.wantInput == False:
                app.wantInput = True
                while True:
                    if app.pdb != '':
                        break
                    try:
                        getUserInput(app,"Please enter a valid 4-digit PDB code\ne.g. 2JXP, 7R6Z, 2FAT, 2P2D, 1SPF\n(Find more on rcsb.org)")
                        break
                    except:
                        pass
            else:
                app.wantInput = False
        elif inButton(event, app.buttonBottomRight):
            goToVisualPage(app)
        elif inButton(event, app.buttonBottomMiddle):
            goToDesignPage(app)
    elif app.isVisualPage:
        if inButton(event, app.buttonSwitch):
            if app.isAtomicModel:
                app.isAtomicModel = False
                app.isSecondaryStruct = True
                app.isRotatingAroundX, app.isRotatingAroundY = False, False
                app.isRotatingAroundZ = True
            else:
                app.isAtomicModel = True
                app.isSecondaryStruct = False
                displaySecondaryStruct(app)
        elif inButton(event, app.buttonTopLeft):
            resetApp(app)
        elif inButton(event, app.buttonX):
            rotateX(app)
        elif inButton(event, app.buttonY):
            rotateY(app)
        elif inButton(event, app.buttonZ):
            rotateZ(app)
        elif inButton(event, app.leftArrowCoords):
            for atom in app.atoms:
                atom.rotateAroundZ(-1)
                atom.coordinate2D = threeDToTwoD(atom.coordinate)
                app.upButtonColor = "SteelBlue1"
                app.downButtonColor = "SteelBlue1"
                app.leftButtonColor = "SteelBlue3"
                app.rightButtonColor = "SteelBlue1"
            app.atoms = sorted(app.atoms, key=cmp_to_key(compare))
        elif inButton(event, app.rightArrowCoords):
            for atom in app.atoms:
                atom.rotateAroundZ(1)
                atom.coordinate2D = threeDToTwoD(atom.coordinate)
                app.upButtonColor = "SteelBlue1"
                app.downButtonColor = "SteelBlue1"
                app.leftButtonColor = "SteelBlue1"
                app.rightButtonColor = "SteelBlue3"
            app.atoms = sorted(app.atoms, key=cmp_to_key(compare))
        elif inButton(event, app.upArrowCoords):
            for atom in app.atoms:
                atom.rotateAroundX(1)
                atom.rotateAroundY(1)
                atom.coordinate2D = threeDToTwoD(atom.coordinate)
                app.upButtonColor = "SteelBlue3"
                app.downButtonColor = "SteelBlue1"
                app.leftButtonColor = "SteelBlue1"
                app.rightButtonColor = "SteelBlue1"
            app.atoms = sorted(app.atoms, key=cmp_to_key(compare))
        elif inButton(event, app.downArrowCoords):
            for atom in app.atoms:
                atom.rotateAroundX(-1)
                atom.rotateAroundY(-1)
                atom.coordinate2D = threeDToTwoD(atom.coordinate)
                app.upButtonColor = "SteelBlue1"
                app.downButtonColor = "SteelBlue3"
                app.leftButtonColor = "SteelBlue1"
                app.rightButtonColor = "SteelBlue1"
            app.atoms = sorted(app.atoms, key=cmp_to_key(compare))

def displaySecondaryStruct(app):
    app.isAtomicModel = True
    app.isSecondaryStruct = False

def mouseDragged(app, event): 
    pass

def userInputToDNA(app, input):
    result = ''
    for i in range(len(input)):
        if input[i] == 'A' or input[i] == 'T' or input[i] == 'G' or input[i] == 'C':
            result += input[i]
    app.input_dna = result

def processUserInput(app):
    pdbcode = app.pdb
    # Sources for retrieving pdb file
    # Source 1: https://www.tutorialspoint.com/biopython/biopython_pdb_module.htm
    # Source 2: https://biopython.org/docs/1.75/api/Bio.PDB.PDBList.html
    PDBList().retrieve_pdb_file(f'{pdbcode}', file_format = 'pdb', pdir = '.',)
    pdb = Path(f'pdb{pdbcode}.ent').read_text()
    coordinates = []
    elements = []
    helices = dict()
    #         {Index: ["{FirstAA}", startPos,  "{LastAA}",  endPos, len]}
    #  e.g.   {1:     ["ASN",          68,         "VAL",    78,    10 ]}
    sheets = dict()
    aminoAcidSeq = dict()
    sheetsCounter = 0
    # Getting amino acid sequence to be converted into RNA & then DNA
    # for user to edit
    aminoAcids = dict()
    AACounter = 1
    for line in pdb.split("\n"):
        if line.startswith("TER"):
            break
        elif line.startswith("ATOM"):
            entries = line.split(" ")
            result = []
            for entry in entries:
                if entry != "":
                    result.append(entry)
            # Amino acid dictionary starts here 
            currAA = result[3]
            if result[1] == '1':
                aminoAcids[AACounter] = currAA
            elif currAA != aminoAcids[AACounter]:
                AACounter += 1
                aminoAcids[AACounter] = currAA
            # Adds to dictionaries of coordinates and elements
            [x, y, z, w] = [5 * float(result[6]), 5 * float(result[7]), 5 * float(result[8]), 1]
            coordinates.append([x,y,z, w])
            elements.append(result[-1])
            # Adds to indexed dictionary of amino acids for modeling secondary 
            # structures
            aaIndex = result[5]
            if aaIndex in aminoAcidSeq:
                [x, y, z, w] = [5 * float(result[6]), 5 * float(result[7]), 5 * float(result[8]), 1]
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
            entries = []
            for entry in line.split(" "):
                if entry != "" and entry != "SHEET":
                    entries.append(entry)
            if len(entries[4]) > 1:
                continue
            elif len(entries) > 17:
                sheetsCounter += 1
                sheets[str(sheetsCounter)] = [entries[3]] + entries[5:7] + [entries[8]]
                sheetsCounter += 1
                sheets[str(sheetsCounter)] = [entries[11]] + [entries[13]] + [entries[15]] + [entries[17]]
            elif len(entries) > 8:
                sheetsCounter += 1
                sheets[str(sheetsCounter)] = [entries[3]] + entries[5:7] + [entries[8]]
    app.aminoAcids = aminoAcids
    app.rna = reverseTranslate(app, AACounter)
    app.dna = reverseTranscribe(app)
    getDNAComplement(app)
    app.aminoAcidSeq = aminoAcidSeq
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
    app.coordinates, app.elems = coordinatesInNumpy, elemsInNumpy
    updateHelices(app, helices)
    app.sheets = sheets
    updateSheets(app, sheets)
    app.atoms = list()
    #            C, H, O, N, P, S
    elemCount = [0, 0, 0, 0, 0, 0]
    for i in range(len(coordinatesInNumpy)):
        app.atoms.append(Protein(app.coordinates[i], app.elems[i]))
        if app.elems[i] == 'C':
            elemCount[0] += 1
        elif app.elems[i] == 'H':
            elemCount[1] += 1
        elif app.elems[i] == 'O':
            elemCount[2] += 1
        elif app.elems[i] == 'N':
            elemCount[3] += 1
        elif app.elems[i] == 'P':
            elemCount[4] += 1
        elif app.elems[i] == 'S':
            elemCount[5] += 1
    totalCount = sum(elemCount)
    app.percent = dict()
    app.percent['Carbon'] = 100 * elemCount[0] / totalCount
    app.percent['Hydrogen'] = 100 * elemCount[1] / totalCount
    app.percent['Oxygen'] = 100 * elemCount[2] / totalCount
    app.percent['Nitrogen'] = 100 * elemCount[3] / totalCount
    app.percent['Phosphorus'] = 100 * elemCount[4] / totalCount
    app.percent['Sulfur'] = 100 * elemCount[5] / totalCount
    for helix in helices:
        startHelixPos, endHelixPos = helices[helix][1], helices[helix][3]
        if startHelixPos in aminoAcidSeq and endHelixPos in aminoAcidSeq:
            startAAPos = aminoAcidSeq[startHelixPos][0]
            endAAPos = aminoAcidSeq[endHelixPos][-1]
            app.structHelix.append(Helix(startAAPos, endAAPos))
    for sheet in sheets:
        startSheetPos, endSheetPos = sheets[sheet][1], sheets[sheet][3]
        if startSheetPos in aminoAcidSeq and endSheetPos in aminoAcidSeq:
            pos = []
            for i in range(int(startSheetPos), int(endSheetPos) + 1):
                if str(i) in aminoAcidSeq:
                    pos.append(aminoAcidSeq[str(i)][0])
            posInTuple = tuple(pos)
            if len(posInTuple) > 0:
                app.structSheet.append(Sheet(posInTuple))
    app.isParsing = False
    app.finishedParsing = True

def updateHelices(app, helices):
    app.helices = helices

def updateSheets(app, sheets):
    app.sheets = sheets

def rotateX(app):
    app.isRotatingX = True
    app.isRotatingY = False
    app.isRotatingZ = False
    app.buttonXColor = "DeepSkyBlue4"
    app.buttonYColor = "DeepSkyBlue2"
    app.buttonZColor = "DeepSkyBlue2"

def rotateY(app):
    app.isRotatingX = False
    app.isRotatingY = True
    app.isRotatingZ = False
    app.buttonXColor = "DeepSkyBlue2"
    app.buttonYColor = "DeepSkyBlue4"
    app.buttonZColor = "DeepSkyBlue2"

def rotateZ(app):
    app.isRotatingX = False
    app.isRotatingY = False
    app.isRotatingZ = True
    app.buttonXColor = "DeepSkyBlue2"
    app.buttonYColor = "DeepSkyBlue2"
    app.buttonZColor = "DeepSkyBlue4"

def goToDesignPage(app):
    app.isIntropage = False
    app.isDesignPage = True
    app.isHomepage = False

def goToPromptWindow(app):
    app.wantInput = True

def goToVisualPage(app):
    app.isHomepage = False
    app.isDesignPage = False
    app.isVisualPage = True
    app.isHelpPage = False
    app.isIntroPage = False 

def redrawAll(app, canvas):
    drawBackground(app, canvas, "mintcream")
    if app.isHomepage:
        drawHomepage(app, canvas)
        drawHelpHint(app, canvas, "snow")
    elif app.isDesignPage:
        drawDesignPage(app, canvas)
        if app.scissors_selected:
            canvas.create_text(app.width/4, 7.5*app.height/8, 
                        text= "Click on a nucleotide to delete DNA bases",
                        fill = "snow", font="Arial 15 bold")
            canvas.create_image(5.2 * app.width / 6, 5.2 * app.height / 6, 
                        image=ImageTk.PhotoImage(app.img_scissors_outline))
        else:
            canvas.create_image(5.2 * app.width / 6, 5.2 * app.height / 6, 
                        image=ImageTk.PhotoImage(app.img_scissors))
        if app.pencil_selected:
            canvas.create_text(app.width/4, 7.5*app.height/8, 
                        text= "Click on a nucleotide to insert DNA",
                        fill = "snow", font="Arial 15 bold")
        else: 
            canvas.create_image(4.8 * app.width / 6, 5.2 * app.height / 6, 
                        image=ImageTk.PhotoImage(app.img_pencil))
        if not app.scissors_selected and not app.pencil_selected:
            canvas.create_text(app.width/4, 7.5*app.height/8, 
                        text= "Select scissors or pencil to edit the DNA sequence",
                        fill = "snow", font="Arial 15 bold")
        drawHelpHint(app, canvas, "RoyalBlue4")
    elif app.isIntroPage:
        drawIntroPage(app, canvas)
        drawHelpHint(app, canvas, "RoyalBlue4")
    elif app.isVisualPage:
        drawVisualPage(app, canvas)
        drawAtomicModel(app, canvas)
        drawSecondaryStruct(app, canvas)
        # Tiny note for buttons
        canvas.create_text(app.buttonX[0] - 10, app.buttonX[1] - 30, 
                        text = "Select axis of rotation", fill = "snow",
                        font = "Arial 15", anchor = "w")
    elif app.isHelpPage:
        drawHelpPage(app,canvas)

def drawCenterButton(app, canvas):
    canvas.create_rectangle(7 * app.width / 8, 7 * app.height / 8, 
                            20 + 7 * app.width / 8, 20 + 7 * app.height /8, 
                            fill = "DeepSkyBlue4")

def drawLeftArrow(app, canvas):
    xOrigin, yOrigin = 6.9 * app.width / 8, 7 * app.height / 8
    canvas.create_polygon(xOrigin, yOrigin,
                          xOrigin - 20, yOrigin,
                          xOrigin - 20, yOrigin - 10,
                          xOrigin - 30, yOrigin + 10,
                          xOrigin - 20, yOrigin + 30,
                          xOrigin - 20, yOrigin + 20,
                          xOrigin, yOrigin + 20, fill = f"{app.leftButtonColor}")

def drawRightArrow(app, canvas):
    xOrigin, yOrigin = 7.1 * app.width / 8 + 20, 7 * app.height / 8
    canvas.create_polygon(xOrigin, yOrigin,
                          xOrigin + 20, yOrigin,
                          xOrigin + 20, yOrigin - 10,
                          xOrigin + 30, yOrigin + 10,
                          xOrigin + 20, yOrigin + 30,
                          xOrigin + 20, yOrigin + 20,
                          xOrigin, yOrigin + 20, fill = f"{app.rightButtonColor}")

def drawDownArrow(app, canvas):
    xOrigin, yOrigin = 7 * app.width / 8, 7.3 * app.height / 8
    canvas.create_polygon(xOrigin, yOrigin,
                          xOrigin, yOrigin + 20,
                          xOrigin - 10, yOrigin + 20,
                          xOrigin + 10, yOrigin + 30,
                          xOrigin + 30, yOrigin + 20,
                          xOrigin + 20, yOrigin + 20,
                          xOrigin + 20, yOrigin, fill = f"{app.downButtonColor}")

def drawUpArrow(app, canvas):
    xOrigin, yOrigin = 7 * app.width / 8, 6.9 * app.height / 8
    canvas.create_polygon(xOrigin, yOrigin,
                          xOrigin, yOrigin - 20,
                          xOrigin - 10, yOrigin - 20,
                          xOrigin + 10, yOrigin - 30,
                          xOrigin + 30, yOrigin - 20,
                          xOrigin + 20, yOrigin - 20,
                          xOrigin + 20, yOrigin, fill = f"{app.upButtonColor}")

def inButton(event, button):
    return (button[0] <= event.x <= button[2] and 
            button[1] <= event.y <= button[3])

def drawHelpPage(app, canvas):
    canvas.create_rectangle(0, 0, app.width, app.height, fill='DarkSeaGreen1')
    canvas.create_text(app.width/2, app.height/2, 
                        text = "press X to return", font = "Arial 20")
    canvas.create_text(app.width / 3, app.height / 5,
                       text="Provide instructions for user here", 
                       font = "Arial 20")

def drawHomepage(app, canvas):
    button1 = app.buttonCenter1
    canvas.create_image(app.width /2 , app.height / 2, 
                        image=ImageTk.PhotoImage(app.image_home))
    canvas.create_text(app.width / 2, 1.8 * app.height / 3,
                    text = "myCrispr", 
                    font = "Helvetica 70", fill="white")
    canvas.create_text(app.width / 2, 2.2 * app.height / 3,
                    text = "CRISPR cas9 enzyme creator at your fingertips", 
                    font = "Arial 16", fill="CadetBlue1")
    if app.center1ButtonMoved:
        drawCenter1Button(app, canvas, f"{app.center1ButtonColor}")    
    else:
        drawCenter1Button(app, canvas, f"{app.center1ButtonColor}")
    canvas.create_text((button1[0] + button1[2])/2, (button1[1] + button1[3])/2,
                       text="Skip intro", fill="Slategray2", 
                       font = "Arial 15 bold")
    canvas.create_text((button1[0] + button1[2])/2, 0.92 * (button1[1] + button1[3])/2,
                       text="press any key to begin", fill="white",
                       font="Arial 20 bold", anchor = "n")

def drawHelpHint(app, canvas, color):
    canvas.create_text(app.width/2, 7.5*app.height/8, text= "press h for help",
                        fill = f"{color}", font="Arial 13")

def drawIntroPage(app, canvas):
    # Source: https://www.nature.com/articles/s41467-018-04252-2
    canvas.create_text(app.width / 10, app.height / 5,
    text="CRISPR-Cas9 is a gene-editing tool and indispensable tool in precision medicine", font = "Arial 20", anchor = "w")
    canvas.create_text(app.width / 10, 1.2*app.height / 5,
    text="and biological research with applications in targeted gene regulation, epigenetic", font = "Arial 20", anchor = "w")
    canvas.create_text(app.width / 10, 1.4*app.height / 5,
    text="modulation, chromatin manipulation, and live cell chromatin imaging.", font = "Arial 20",  anchor = "w")
    canvas.create_text(app.width / 10, 1.7*app.height / 5,
    text="Designing your own CRISPR enzyme is no rocket science.", font = "Arial 20",  anchor = "w")
    canvas.create_text(app.width / 10, 1.9*app.height / 5,
    text="Find your favorite genetic sequence and create your own CRSIPR-Cas9 enzyme using myCrispr.", font = "Arial 20",  anchor = "w")
    drawReturnButton(app,canvas)
    drawBottomLeftButton(app, canvas, "steel blue")
    canvas.create_image(app.width / 3.3, app.height / 1.7,
                        image=ImageTk.PhotoImage(app.image_design))
    buttonLeft = app.buttonBottomLeft
    canvas.create_text((buttonLeft[0] + buttonLeft[2])/2, (buttonLeft[1] + buttonLeft[3])/2,
                       text = "Choose Your Protein", fill = "snow", 
                       font = "Arial 15 bold")
    buttonRight = app.buttonBottomRight
    if app.pdb != '':
        canvas.create_image(2 * app.width / 2.85, app.height / 1.7,
                            image=ImageTk.PhotoImage(app.protein_model))
        drawBottomRightButton(app, canvas, "pale violet red")
        canvas.create_text((buttonRight[0] + buttonRight[2])/2, (buttonRight[1]+
                            buttonRight[3])/2, text = "Visualize Protein", 
                            fill = "snow", font = "Arial 15 bold")
        drawBottomMiddleButton(app, canvas, "RosyBrown3")
        buttonMiddle = app.buttonBottomMiddle
        canvas.create_text((buttonMiddle[0] + buttonMiddle[2])/2, (buttonMiddle[1]+
                            buttonMiddle[3])/2, text = "Edit Protein", 
                            fill = "grey1", font = "Arial 15 bold")

def drawDesignPage(app, canvas):
    drawBackground(app, canvas, "#001a33")
    drawReturnButton(app,canvas)
    # Draw a triangle indicating the current base
    canvas.create_polygon(85, app.y2 + 30, 75, app.y2 + 50, 95, app.y2 + 50,
                          fill="yellow")
    currNum = getIndex(app, 85)
    totalNum = len(app.dna)
    canvas.create_text(85, app.y2 + 60, text = f"Current Base: {currNum}/{totalNum}",
                       fill = "snow", font = "Arial 15 bold")
    if currNum < totalNum * .1:
        canvas.create_text(120, app.y2 + 80, text = "Press right arrow key to move",
                       fill = "snow", font = "Arial 15 bold")
    elif currNum > totalNum * .1 and currNum < totalNum * .9:
        canvas.create_text(140, app.y2 + 80, text = "Press left or right arrow key to move",
                       fill = "snow", font = "Arial 15 bold")
    else:
        canvas.create_text(120, app.y2 + 80, text = "Press left arrow key to move",
                       fill = "snow", font = "Arial 15 bold")
    # Draw neon yellow highlight on dna double helix
    ratio = getIndex(app, 0) / len(app.dna)
    tPos = ratio * 1000
    lower = 100 * np.sin(2 * np.pi * (tPos + 240) / 255) + app.height/3.5
    upper = 100 * np.cos(2 * np.pi * (tPos + 15) / 255) + app.height/3.5
    canvas.create_line(tPos, lower, tPos, upper, fill = "yellow", width = 10)
    seq = app.dna
    spacing = 1000 / len(seq)
    # Draw double helix
    for i in range(len(seq)):
        t = i * spacing 
        x = 100 * np.sin(2 * np.pi * (t + 240) / 255) + app.height/3.5
        y = 100 * np.cos(2 * np.pi * (t + 15) / 255) + app.height/3.5
        # if x < y:
        #     canvas.create_line(t, x, t, y)
        # else: 
        #     canvas.create_line(t, y, t, x)
    for t in app.timesteps:
        # RGB values range from 0 to 255
        #                    dark -> light
        # Timesteps range from 0 to 1000
        r = 12
        # First strand of sugar-phosphate backbone
        x = 100 * np.sin(2 * np.pi * (t + 240) / 255) + app.height/3.5
        canvas.create_oval(t-r, x-r, t+r, x+r, fill = f"#{toHex(toRGB(t))}",
                                outline = f"#{toHex(toRGB(t))}")
        # Second strand of sugar-phosphate backbone
        y = 100 * np.cos(2 * np.pi * (t + 15) / 255) + app.height/3.5
        canvas.create_oval(t-r, y-r, t+r, y+r, fill = f"#{toHex(toRGB(t))}", 
                            outline=f"#{toHex(toRGB(t))}")
    numOfBlocks = len(app.dna)
    y0, y1, y2 = app.y0, app.y1, app.y2
    # Draw lines that show the scale
    canvas.create_line(tPos, max(lower, upper), 0, y0, fill = "light sky blue", width = 5)
    canvas.create_line(tPos, max(lower, upper), 1000, y0, fill = "light sky blue", width = 5)
    scale = app.blockZoomScale
    spacing = app.spacing
    currSpacing = 0
    # Draw DNA nucleotides
    for i in range(numOfBlocks):
        currSpacing = i * spacing
        x0, x1 = app.scrollX + currSpacing + i * scale, app.scrollX + currSpacing + (i+1) * scale
        if app.dna[i] == 'A':
            blockColor = "light salmon"
            compColor = "lemon chiffon"
            textColor = "grey3"
        elif app.dna[i] == 'T':
            blockColor = "lemon chiffon"
            compColor = "light salmon"
            textColor = "grey3"
        elif app.dna[i] == 'G':
            blockColor = "LightBlue1"
            compColor = "DarkSeaGreen1"
            textColor = "grey3"
        elif app.dna[i] == 'C':
            blockColor = "DarkSeaGreen1"
            compColor = "LightBlue1"
            textColor = "grey3"
        canvas.create_rectangle(x0, y0, x1, y1, fill=f"{blockColor}")
        canvas.create_text((x0+x1)/2, (y0+y1)/2, text=f"{app.dna[i]}",
                            fill=f"{textColor}")
        canvas.create_rectangle(x0, y1, x1, y2, fill=f"{compColor}")
        canvas.create_text((x0+x1)/2, (y1+y2)/2, text=f"{app.DNAComplement[i]}",
                            fill=f"{textColor}")
    # Draw highlighter box for selected nucleotide
    if app.isDrawingHighlighterBox:
        i = app.currIndex
        x0 = app.scrollX + i * (scale + spacing)
        x1 = app.scrollX + (i+1) * (scale + spacing)
        canvas.create_rectangle(x0 - 3, y0 - 3, x1 - 2, y2 + 3, outline="yellow", width = 3)

def getIndex(app, x0):
    index = (x0 - app.scrollX) // (app.spacing + app.blockZoomScale)
    return int(index)

def toRGB(t):
    r = int(abs(150 * np.sin(2 * np.pi * (t-200) / 255)))
    g = int(0.8 * r)
    b = 255
    return (r, g, b)

# Source of the toHex() function:
# https://www.codespeedy.com/convert-rgb-to-hex-color-code-in-python/
def toHex(tuple):
    return "%02x%02x%02x" % tuple

def drawVisualPage(app, canvas):
    drawBackground(app, canvas, "grey8")
    drawReturnButton(app, canvas)
    drawAxes(app, canvas)
    drawButtons(app, canvas)

def drawButtons(app, canvas):
    # Button X
    canvas.create_rectangle(app.buttonX[0], app.buttonX[1], app.buttonX[2],
                            app.buttonX[3],
                            fill = f"{app.buttonXColor}")
    canvas.create_text((app.buttonX[0] + app.buttonX[2])/2,
                        (app.buttonX[1] + app.buttonX[3])/2,
                        text = "X", font = "Arial 15 bold", fill = "snow")
    # Button Y
    canvas.create_rectangle(app.buttonY[0], app.buttonY[1], app.buttonY[2],
                            app.buttonY[3],
                            fill = f"{app.buttonYColor}")
    canvas.create_text((app.buttonY[0] + app.buttonY[2])/2,
                        (app.buttonY[1] + app.buttonY[3])/2,
                        text = "Y", font = "Arial 15 bold", fill = "snow")
    # Button Z
    canvas.create_rectangle(app.buttonZ[0], app.buttonZ[1], app.buttonZ[2],
                            app.buttonZ[3],
                            fill = f"{app.buttonZColor}")
    canvas.create_text((app.buttonZ[0] + app.buttonZ[2])/2,
                        (app.buttonZ[1] + app.buttonZ[3])/2,
                        text = "Z", font = "Arial 15 bold", fill = "snow")

def drawAtomicModel(app, canvas):
    if app.isAtomicModel:
        for atom in app.atoms:
            coordinate = (atom.coordinate2D[0] + app.width / 2, 
                        atom.coordinate2D[1] + app.width / 3)
            if atom.atom == 'C':
                drawCarbon(app, canvas, coordinate, atom)
            elif atom.atom == 'H':
                drawHydrogen(app, canvas, coordinate, atom)
            elif atom.atom == 'O':
                drawOxygen(app, canvas, coordinate, atom)
            elif atom.atom == 'N':
                drawNitrogen(app, canvas, coordinate, atom)
            elif atom.atom == 'P':
                drawPhosphorus(app, canvas, coordinate, atom)
            elif atom.atom == 'S':
                drawSulfur(app, canvas, coordinate, atom)
        # Displays chemical composition
        counter = 0
        for element in app.percent:
            canvas.create_text(3.5 * app.width / 5, (1 + 0.3 * counter) * app.height / 10, 
                        text = f"{element}: {app.percent[element]}%", anchor = "w", 
                        fill = "mint cream", font = "Arial 15 bold")
            counter += 1
        # Draw arrows for the user to rotate the protein
        drawCenterButton(app, canvas)
        drawLeftArrow(app, canvas)
        drawRightArrow(app, canvas)
        drawUpArrow(app, canvas)
        drawDownArrow(app, canvas)
        drawSwitchButton(app, canvas)

def drawSecondaryStruct(app, canvas):
    helixStructs = app.structHelix
    sheetStructs = app.structSheet
    if app.isSecondaryStruct:
        for i in range(len(sheetStructs)):
            struct = sheetStructs[i]
            posTuple = struct.posIn2D
            color = app.sheetColors[i % 9]
            drawSheet(app, canvas, posTuple, color)
        for struct in helixStructs:
            startPos, endPos = struct.start, struct.end
            drawHelix(app, canvas, startPos, endPos)
    drawCenterButton(app, canvas)
    drawLeftArrow(app, canvas)
    drawRightArrow(app, canvas)
    drawUpArrow(app, canvas)
    drawDownArrow(app, canvas)
    drawSwitchButton(app, canvas)

def drawHelix(app, canvas, start, end):
    startPos, endPos = threeDToTwoD(start), threeDToTwoD(end)
    x0, x1 = startPos[0] + app.width / 3, endPos[0] + app.width / 3
    y0, y1 = startPos[1] + app.height / 4 , endPos[1] + app.height / 4
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

def drawSheet(app, canvas, posTuple, color):
    # startPos, endPos = threeDToTwoD(start), threeDToTwoD(end)
    # x0, x1 = startPos[0] + app.width / 3, endPos[0] + app.width / 3
    # y0, y1 = startPos[1] + app.height / 4 , endPos[1] + app.height / 4
    # # The original image is 1998 pixels wide (we only care about the length)
    # dx, dy = x1 - x0, y1 - y0
    if len(posTuple) >= 2:
        canvas.create_line(posTuple, width = 10, fill=f"{color}")
        # canvas.create_line(posTuple[0][0], posTuple[-1][0],
        #                    posTuple[0][1], posTuple[-1][1], 
        #                    width = 10, fill = "yellow")
    # distance = (dx ** 2 + dy ** 2) ** 0.5
    # ratio = distance / 1998
    # image = app.image_sheet
    # radOfRot = math.atan(dy/dx)
    # degreeOfRot = math.degrees(radOfRot)
    # imageOfSheet = image.rotate(degreeOfRot)
    # # The ratio is further scaled since the original ratio would render 
    # # the helices too large
    # image_sheet_scaled = app.scaleImage(imageOfSheet, ratio)
    # # Anchor the center of the image to the midpoint between the start 
    # # and end pos
    # center = ((x0 + x1) / 2, (y0 + y1) / 2)
    # canvas.create_image(center[0], center[1], 
    #                     image=ImageTk.PhotoImage(image_sheet_scaled))

def drawBottomLeftButton(app, canvas, color):
    button = app.buttonBottomLeft
    canvas.create_rectangle(button[0], button[1], button[2], button[3],
                            fill=f"{color}")

def drawSwitchButton(app, canvas):
    button = app.buttonSwitch
    canvas.create_rectangle(button[0], button[1], button[2], button[3],
                            fill=f"{app.buttonSwitchColor}")
    if app.isAtomicModel:
        canvas.create_text((button[0]+button[2])/2, (button[1]+button[3])/2,
                            text="Secondary\nStructures", font = "Arial 10")
    elif app.isSecondaryStruct:
        canvas.create_text((button[0]+button[2])/2, (button[1]+button[3])/2,
                            text="Atomic View", font = "Arial 10")

def drawBottomRightButton(app, canvas, color):
    button = app.buttonBottomRight
    canvas.create_rectangle(button[0], button[1], button[2], button[3],
                            fill=f"{color}")

def drawBottomMiddleButton(app, canvas, color):
    button = app.buttonBottomMiddle
    canvas.create_rectangle(button[0], button[1], button[2], button[3],
                            fill=f"{color}")

def drawCenter1Button(app, canvas, color):
    button1 = app.buttonCenter1
    canvas.create_rectangle(button1[0], button1[1], button1[2], button1[3],
                            fill = f"{color}")

def drawReturnButton(app,canvas):
    buttonReturn = app.buttonTopLeft
    canvas.create_rectangle(buttonReturn[0], buttonReturn[1], buttonReturn[2], 
                            buttonReturn[3], fill = "RoyalBlue3")
    canvas.create_text((buttonReturn[0]+buttonReturn[2])/2, 
                        (buttonReturn[1]+buttonReturn[3])/2,
                        text = "Return", font = "Arial 20", fill = "snow")

def drawBackground(app, canvas, color):
    canvas.create_rectangle(0, 0, app.width, app.height, fill = f"{color}")

def drawAxes(app, canvas):
    # X axis
    canvas.create_line(app.width / 5, 4 * app.height / 5, app.width / 2, app.height / 2)
    # Y axis
    canvas.create_line(app.width / 2, app.height / 2, 4 * app.width / 5, 4 * app.height / 5)
    # Z axis
    canvas.create_line(app.width / 2, app.height / 5, app.width / 2, app.height / 2)

# Atomic radii of the elements
# C =  77, H = 37, O = 73, N = 70, P = 110, S = 103
def drawCarbon(app, canvas, coordinate, atom):
    coords = atom.coordinate[:3]
    myCoords = app.viewerCoords
    rgbScalar = int(getNorm(myCoords[0] - coords[0], myCoords[1] - coords[1], myCoords[2] - coords[2])/5)
    r = 0.9 * 77/73 + 0.05 * app.level
    (x, y) = (coordinate[0], coordinate[1])
    canvas.create_oval(x-r, y-r, x+r, y+r, fill = f'#{toHex((rgbScalar, rgbScalar, rgbScalar))}')

def drawHydrogen(app, canvas, coordinate, atom):
    r = 37/73 + 0.05 * app.level
    (x, y) = (coordinate[0], coordinate[1])
    canvas.create_oval(x-r, y-r, x+r, y+r, fill = 'floral white')

def drawOxygen(app, canvas, coordinate, atom):
    r = 73/73 + 0.05 * app.level
    (x, y) = (coordinate[0], coordinate[1])
    canvas.create_oval(x-r, y-r, x+r, y+r, fill = 'IndianRed3')

def drawNitrogen(app, canvas, coordinate, atom):
    r = 70/73 + 0.05 * app.level
    (x, y) = (coordinate[0], coordinate[1])
    canvas.create_oval(x-r, y-r, x+r, y+r, fill = 'RoyalBlue3')

def drawPhosphorus(app, canvas, coordinate, atom):
    r = 110/73 + 0.05 * app.level
    (x, y) = (coordinate[0], coordinate[1])
    canvas.create_oval(x-r, y-r, x+r, y+r, fill = 'gold', outline = "gold")

def drawSulfur(app, canvas, coordinate, atom):
    r = 1.5 * 103/73 + 0.05 * app.level
    (x, y) = (coordinate[0], coordinate[1])
    canvas.create_oval(x-r, y-r, x+r, y+r, fill = 'LightPink1')

def getNorm(x, y, z):
    return (x ** 2 + y ** 2 + z ** 2) ** 0.5

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

class Sheet(object):
    def __init__(self, posInTuple):
        self.pos = posInTuple
        result = []
        for tempPos in posInTuple:
            realPos = tuple(threeDToTwoD(tempPos))
            # app.width, app.height = 1000, 800
            adjustedPos = (realPos[0] + 1000 / 3, realPos[1] + 800 / 4)
            result.append(adjustedPos)
        self.posIn2D = tuple(result)

    def scaleUp(self):
        entries = [[1.05, 0, 0, 0], 
                   [0, 1.05, 0, 0],
                   [0, 0, 1.05, 0],
                   [0, 0, 0, 1]]
        scaleUpMatrix = np.array(entries)
        result = []
        for tempPos in self.pos:
            newPos = np.matmul(scaleUpMatrix, tempPos)
            result.append(tuple(newPos))
        self.pos = tuple(result)

    def scaleDown(self):
        entries = [[0.95, 0, 0, 0], 
                   [0, 0.95, 0, 0],
                   [0, 0, 0.95, 0],
                   [0, 0, 0, 1]]
        scaleDownMatrix = np.array(entries)
        result = []
        for tempPos in self.pos:
            newPos = np.matmul(scaleDownMatrix, tempPos)
            result.append(tuple(newPos))
        self.pos = tuple(result)

    def rotateAroundX(self, sign):
        entries = [[math.cos(sign * math.pi / 12), -math.sin(sign * math.pi / 12), 0, 0],
                   [math.sin(sign * math.pi / 12), math.cos(sign * math.pi / 12), 0, 0],
                   [0, 0, 1, 0],
                   [0, 0, 0, 1]]
        rotateMatrix = np.array(entries)
        result = []
        for tempPos in self.pos:
            newPos = np.matmul(rotateMatrix, tempPos)
            result.append(newPos)
        self.pos = tuple(result)
        result = []
        for temp2DPos in self.pos:
            new2DPos = tuple(threeDToTwoD(temp2DPos))
            shifted2DPos = (new2DPos[0] + 1000 / 3, new2DPos[1] + 800 / 4)
            result.append(tuple(shifted2DPos))
        self.posIn2D = tuple(result)

    def rotateAroundY(self, sign):
        entries = [[1, 0, 0, 0],
                   [0, math.cos(sign * math.pi / 12), math.sin(sign * math.pi / 12), 0],
                   [0, -math.sin(sign * math.pi / 12), math.cos(sign * math.pi / 12), 0],
                   [0, 0, 0, 1]]
        rotateMatrix = np.array(entries)
        result = []
        for tempPos in self.pos:
            newPos = np.matmul(rotateMatrix, tempPos)
            result.append(newPos)
        self.pos = tuple(result)
        result = []
        for temp2DPos in self.pos:
            new2DPos = tuple(threeDToTwoD(temp2DPos))
            shifted2DPos = (new2DPos[0] + 1000 / 3, new2DPos[1] + 800 / 4)
            result.append(tuple(shifted2DPos))
        self.posIn2D = tuple(result)

    def rotateAroundZ(self, sign):
        entries = [[math.cos(sign * math.pi / 12), 0, math.sin(sign * math.pi / 12), 0],
                   [0, 1, 0, 0],
                   [-math.sin(sign * math.pi / 12), 0, math.cos(sign * math.pi / 12), 0],
                   [0, 0, 0, 1]]
        rotateMatrix = np.array(entries)
        result = []
        for tempPos in self.pos:
            newPos = np.matmul(rotateMatrix, tempPos)
            result.append(newPos)
        self.pos = tuple(result)
        result = []
        for temp2DPos in self.pos:
            new2DPos = tuple(threeDToTwoD(temp2DPos))
            shifted2DPos = (new2DPos[0] + 1000 / 3, new2DPos[1] + 800 / 4)
            result.append(tuple(shifted2DPos))
        self.posIn2D = tuple(result)

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

def translate(app, newRNA):
    # Translate (RNA -> Protein)
    newProtein = dict()
    nucleotideCount = len(newRNA) // 3 + 1
    for counter in range(nucleotideCount + 1):
        baseStart, baseEnd = counter * 3, counter * 3 + 2
        currNucleotide = newRNA[baseStart:baseEnd+1]
        for AA in app.amino_acids:
            if currNucleotide in app.amino_acids[AA]:
                newProtein[counter + 1] = AA
                break
        if AA == "STOP":
            break
    print(newProtein)
    if len(newProtein) > nucleotideCount:
        app.nonsense_mutation = True
    return newProtein

def transcribe(DNA):
    # Transcribe (DNA -> RNA)
    newRNA = ''
    for baseIndex in range(len(DNA)):
        if DNA[baseIndex] == "T":
            newRNA += "U"
        else:
            newRNA += DNA[baseIndex]
    return newRNA

def reverseTranscribe(app):
    RNA = app.rna
    # RNA -> DNA
    DNA = ''
    for baseIndex in range(len(RNA)):
        if RNA[baseIndex] == "U":
            DNA += "T"
        else:
            DNA += RNA[baseIndex]
    return DNA

def reverseTranslate(app, AACounter):
    # Protein -> RNA
    RNA = ''
    for aminoAcidIndex in range(1, AACounter + 1):
        currAA = app.aminoAcids[aminoAcidIndex]
        codon = random.choice(tuple(app.amino_acids[currAA]))
        RNA += codon
    return RNA

def getDNAComplement(app):
    seq = app.dna
    result = ''
    for i in range(len(seq)):
        if seq[i] == 'A':
            complement = 'T'
        elif seq[i] == 'T':
            complement = 'A'
        elif seq[i] == 'G':
            complement = 'C'
        elif seq[i] == 'C':
            complement = 'G'
        result += complement
    app.DNAComplement = result

runApp(width = 1000, height = 800)