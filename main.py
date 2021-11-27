# Name: Claire Chen
# Andrew ID: ccz
from cmu_112_graphics import *
from seqAnalysis import *
import numpy as np
from projections import *
import math
import numpy as np
import parser
from functools import cmp_to_key
from PIL import Image
from Bio.PDB import *
from pathlib import Path

def appStarted(app):
    app.pdb = ''
    app.isParsing = False
    app.finishedParsing = False
    app.inputs = []
    app.atoms = list()
    app.level = 50
    app.viewerCoords = (500, 500, 500)
    app.coordinates, app.elems = [], []
    app.aminoAcidSeq, app.helices, app.sheets = [], [], []
    url_home = 'https://mma.prnewswire.com/media/1424831/CRISPeR_Gene_Editing.jpg?w=1200'
    app.image_home = app.loadImage(url_home)
    app.image_home_scaled = app.scaleImage(app.image_home, 2/3)
    url_design_scissors = "https://www.sciencenewsforstudents.org/wp-content/uploads/2019/11/080819_ti_crisprsicklecell_feat-1028x579-1028x579.jpg"
    image_design_scissors = app.loadImage(url_design_scissors)
    mouse = app.scaleImage(image_design_scissors, 1/3)
    app.image_mouse = mouse.crop((200, 200, 400, 579))
    app.buttonCenter1 = (2 * app.width / 5, 4.15 * app.height / 5,
                         3 * app.width / 5,  4.4 * app.height / 5)
    app.buttonTopLeft = (app.width / 12, app.height / 12, 
                         app.width / 4, app.height / 8)
    app.buttonBottomLeft= (app.width/5, 3 * app.height / 4, 2 * app.width /5, 
                            3.2 * app.height/4)
    app.buttonBottomRight = (3 * app.width/5, 3 * app.height / 4, 4 * app.width/5, 
                            3.2 * app.height/4)
    app.buttonSettings = (app.width / 10, app.height / 8, 
                          1.1 * app.width / 10, 1.08 * app.height / 8,)
    app.buttonInput = (app.width / 6, 7 * app.height / 8,
                            1.2 * app.width / 6, 7.5 * app.height / 8)
    app.dna_bases = {"adenine", "guanine", "cytosine", "thymine"}
    app.rna_bases = {"adenine", "guanine", "cytosine", "uracil"}
    app.amino_acids = {"ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", 
                       "HIS", "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER",
                       "THR", "TRP", "TYR", "VAL"}
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
    resetApp(app)

def timerFired(app):
    if (app.isVisualPage and app.isAtomicModel and not app.isRotatingX and 
        not app.isRotatingY and not app.isRotatingZ):
        for atom in app.atoms:
            atom.rotateAroundZ(-.3)
            atom.coordinate2D = threeDToTwoD(atom.coordinate)
        # Source: https://stackoverflow.com/questions/5213033/sort-a-list-of-lists-with-a-custom-compare-function
        app.atoms = sorted(app.atoms, key=cmp_to_key(compare))

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
    # Homepage
    app.isHomepage = True
    app.isDesignPage = False
    app.isHelpPage = False
    app.isIntroPage = False 
    app.drawVisualPage = False
    app.drawCustomSeqPage = False
    app.isVisualPage = False

def keyPressed(app, event):
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
                        struct.start2D = threeDToTwoD(struct.start)
                        struct.end2D = threeDToTwoD(struct.end)
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
                        struct.start2D = threeDToTwoD(struct.start)
                        struct.end2D = threeDToTwoD(struct.end)
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
                        struct.start2D = threeDToTwoD(struct.start)
                        struct.end2D = threeDToTwoD(struct.end)
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
                        struct.start2D = threeDToTwoD(struct.start)
                        struct.end2D = threeDToTwoD(struct.end)
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
                        struct.start2D = threeDToTwoD(struct.start)
                        struct.end2D = threeDToTwoD(struct.end)
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
                        struct.start2D = threeDToTwoD(struct.start)
                        struct.end2D = threeDToTwoD(struct.end)
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
                        struct.start2D = threeDToTwoD(struct.start)
                        struct.end2D = threeDToTwoD(struct.end)
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
                    struct.start2D = threeDToTwoD(struct.start)
                    struct.end2D = threeDToTwoD(struct.end)
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
                    struct.rotateAroundY(1)
                    struct.start2D = threeDToTwoD(struct.start)
                    struct.end2D = threeDToTwoD(struct.end)
    elif event.key == "h":
        app.isHelpPage = True
    elif app.isHelpPage and event.key == "x":
        app.isHelpPage = False
    elif app.isHomepage:
        app.isHomepage = False
        app.isIntroPage = True

def getUserInput(app, prompt):
    if app.wantInput:
        # if app.getUserInput(prompt) != None:
        app.inputs.append(app.getUserInput(prompt).lower())
        if len(app.inputs) > 0 and len(app.inputs[-1]) and app.inputs[-1].isalnum():
            app.pdb = app.inputs[-1]
            processUserInput(app)
            app.isParsing = True
            app.wantInput = False
        # else:
        #     return False
    print(app.inputs)

def mousePressed(app, event):
    print(f"app.isAtomicModel is {app.isAtomicModel}")
    print(f"app.isSecondaryStruct is {app.isSecondaryStruct}")
    if app.isHomepage:
        if inButton(event, app.buttonCenter1):
            app.isHomepage = False
            app.isDesignPage = True
    elif app.isDesignPage:
        if inButton(event, app.buttonTopLeft):
            resetApp(app)
    elif app.isIntroPage:
        if inButton(event, app.buttonTopLeft):
            resetApp(app)
        elif inButton(event, app.buttonBottomLeft):
            if app.wantInput == False:
                app.wantInput = True
                getUserInput(app,"Please enter a valid 4-digit PDB code\ne.g. 1S5L, 2FAT, 3ICB")
            else:
                app.wantInput = False
            # goToDesignPage(app)
        elif inButton(event, app.buttonBottomRight):
            goToVisualPage(app)
    elif app.isVisualPage:
        if inButton(event, app.buttonSwitch):
            if app.isAtomicModel:
                app.isAtomicModel = False
                app.isSecondaryStruct = True
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
    print(event)

def displaySecondaryStruct(app):
    app.isAtomicModel = True
    app.isSecondaryStruct = False

def mouseDragged(app, event): 
    pass

def processUserInput(app):
    pdbcode = app.pdb
    # Source 1: https://www.tutorialspoint.com/biopython/biopython_pdb_module.htm
    # Source 2: https://biopython.org/docs/1.75/api/Bio.PDB.PDBList.html
    try:
        PDBList().retrieve_pdb_file(f'{pdbcode}', file_format = 'pdb', pdir = '.',)
        pdb = Path(f'pdb{pdbcode}.ent').read_text()
        print("file downloaded!")
    except: 
        pass

    coordinates = []
    elements = []
    helices = dict()
    #         {Index: ["{FirstAA}", startPos,  "{LastAA}",  endPos, len]}
    #  e.g.   {1:     ["ASN",          68,         "VAL",    78,    10 ]}
    sheets = dict()
    aminoAcidSeq = dict()
    sheetsCounter = 0
    for line in pdb.split("\n"):
        if line.startswith("TER"):
            break
        elif line.startswith("ATOM"):
            entries = line.split(" ")
            result = []
            for entry in entries:
                if entry != "":
                    result.append(entry)
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
    print(aminoAcidSeq)
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
    aminoAcidSeq = app.aminoAcidSeq
    print("helices are", helices)
    app.helices = helices
    app.sheets = sheets
    print("sheets are ", sheets)
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
    print("helices are", helices)
    for helix in helices:
        startHelixPos, endHelixPos = helices[helix][1], helices[helix][3]
        if startHelixPos in aminoAcidSeq and endHelixPos in aminoAcidSeq:
            startAAPos = aminoAcidSeq[startHelixPos][0]
            endAAPos = aminoAcidSeq[endHelixPos][-1]
            # startPos = threeDToTwoD(startAAPos)
            # endPos = threeDToTwoD(endAAPos)
    #         app.structHelix.append((startPos, endPos))
            app.structHelix.append(Helix(startAAPos, endAAPos))
            print("appending helix to app.structHelix...")
    for sheet in sheets:
        startSheetPos, endSheetPos = sheets[sheet][1], sheets[sheet][3]
        if startSheetPos in aminoAcidSeq and endSheetPos in aminoAcidSeq:
            startAAPos = aminoAcidSeq[startSheetPos][0]
            endAAPos = aminoAcidSeq[endSheetPos][-1]
            # startPos = threeDToTwoD(startAAPos)
            # endPos = threeDToTwoD(endAAPos)
            # app.structSheet.append((startPos, endPos))
            app.structSheet.append(Sheet(startAAPos, endAAPos))
    app.isParsing = False
    app.finishedParsing = True

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
        drawHelpHint(app, canvas, "RoyalBlue4")
    elif app.isIntroPage:
        drawIntroPage(app, canvas)
        drawHelpHint(app, canvas, "RoyalBlue4")
    elif app.isVisualPage:
        drawVisualPage(app, canvas)
        drawAtomicModel(app, canvas)
        if app.isSecondaryStruct:
            drawSecondaryStruct(app, canvas)
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

def drawDesignPage(app, canvas):
    drawBackground(app, canvas, "#001a33")
    drawReturnButton(app,canvas)
    seq = sample_dna_info['Sequence']
    spacing = 1000 / len(seq)
    canvas.create_line(0, 400, 1000, 400)
    for i in range(len(seq)):
        t = i * spacing 
        x = 100 * np.sin(2 * np.pi * (t + 240) / 255) + app.height/2
        y = 100 * np.cos(2 * np.pi * (t + 15) / 255) + app.height / 4
        if x < y:
            canvas.create_line(t, x, t, y)
        else: 
            canvas.create_line(t, y, t, x)
    for t in app.timesteps:
        # RGB values range from 0 to 255
        #                    dark -> light
        # Timesteps range from 0 to 1000
        r = 12
        # First strand of sugar-phosphate backbone
        x = 100 * np.sin(2 * np.pi * (t + 240) / 255) + app.height/2
        canvas.create_oval(t-r, x-r, t+r, x+r, fill = f"#{toHex(toRGB(t))}",
                                outline = f"#{toHex(toRGB(t))}")
        # Second strand of sugar-phosphate backbone
        y = 100 * np.cos(2 * np.pi * (t + 15) / 255) + app.height/2
        canvas.create_oval(t-r, y-r, t+r, y+r, fill = f"#{toHex(toRGB(t))}", 
                            outline=f"#{toHex(toRGB(t))}")
    for i in range(len(seq)):
        t = i * spacing 
        x = 100 * np.sin(2 * np.pi * (t + 240) / 255) + app.height/2
        y = 100 * np.cos(2 * np.pi * (t + 15) / 255) + app.height/2
        # Prints original strand
        canvas.create_text(t, (x+400)/2, text = f"{seq[i]}", fill = "mintcream",
                        font = "Arial 15 bold")
        # Prints complementary strand
        if seq[i] == 'A':
            complement = 'T'
        elif seq[i] == 'T':
            complement = 'A'
        elif seq[i] == 'G':
            complement = 'C'
        elif seq[i] == 'C':
            complement = 'G'
        canvas.create_text(t, (y+400)/2, text = f"{complement}", 
                          fill = "mintcream", font = "Arial 15 bold")

def toRGB(t):
    r = int(abs(150 * np.sin(2 * np.pi * (t-200) / 255)))
    g = int(0.8 * r)
    b = 255
    return (r, g, b)

def toHex(tuple):
    return "%02x%02x%02x" % tuple

def drawVisualPage(app, canvas):
    drawBackground(app, canvas, "grey8")
    drawReturnButton(app, canvas)
    drawAxes(app, canvas)

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
        canvas.create_text(app.width / 6, 5 * app.height / 6, text = f"app level = {app.level}")
        # Displays chemical composition
        counter = 0
        for element in app.percent:
            canvas.create_text(3.5 * app.width / 5, (1 + 0.3 * counter) * app.height / 10, 
                        text = f"{element}: {app.percent[element]}%", anchor = "w", 
                        fill = "mint cream", font = "Arial 15 bold")
            counter += 1
        # Tiny note for buttons
        canvas.create_text(app.buttonX[0] - 10, app.buttonX[1] - 30, 
                        text = "Select axis of rotation", fill = "snow",
                        font = "Arial 15", anchor = "w")
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
        # Draw arrows for the user to rotate the protein
        drawCenterButton(app, canvas)
        drawLeftArrow(app, canvas)
        drawRightArrow(app, canvas)
        drawUpArrow(app, canvas)
        drawDownArrow(app, canvas)
        drawSwitchButton(app, canvas)

def drawSecondaryStruct(app, canvas):
    for struct in app.structHelix:
        startPos, endPos = struct.start, struct.end
        print(startPos, endPos)
        # startPos, endPos = struct[0], struct[1]
        drawHelix(app, canvas, startPos, endPos)
    for struct in app.structSheet:
        startPos, endPos = struct.start, struct.end
        # startPos, endPos = struct[0], struct[1]
        drawSheet(app, canvas, startPos, endPos)
    drawCenterButton(app, canvas)
    drawLeftArrow(app, canvas)
    drawRightArrow(app, canvas)
    drawUpArrow(app, canvas)
    drawDownArrow(app, canvas)
    drawSwitchButton(app, canvas)

def drawHelix(app, canvas, startPos, endPos):
    x0, x1 = startPos[0] + app.width / 1.8, endPos[0] + app.width / 1.8
    y0, y1 = startPos[1], endPos[1]
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
    y0, y1 = startPos[1], endPos[1]
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

def drawCustomSeqPage(app, canvas):
    drawReturnButton(app, canvas)
    dna_dict = sample_dna_info
    counter = 0
    for entry in sample_dna_info:
        # Displays dna sequence
        canvas.create_text(2 * app.width / 3, (1 + 0.1 * counter) * app.height /3, 
                        text = f"{entry}", anchor = "w", font = "Arial 15 bold")
        counter += 1
        # Displays protein sequence
        canvas.create_text(2 * app.width / 3, (1 + 0.1 * counter) * app.height /3, 
                        text = dna_dict[f"{entry}"], anchor = "w", font = "Arial 15")
        counter += 1

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

runApp(width = 1000, height = 800)