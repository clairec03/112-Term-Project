# Name: Claire Chen
# Andrew ID: ccz
from cmu_112_graphics import *
from seqAnalysis import *
import numpy as np
from projections import *
import math
import numpy as np
import parser

coordinatesInNumpy = parser.coordinatesInNumpy
elemsInNumpy = parser.elemsInNumpy

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

    def rotateAroundZ(self, sign):
        entries = [[math.cos(sign * math.pi / 12), 0, math.sin(sign * math.pi / 12), 0],
                    [0, 1, 0, 0],
                    [-math.sin(sign * math.pi / 12), 0, math.cos(sign * math.pi / 12), 0],
                    [0, 0, 0, 1]]
        rotateMatrix = np.array(entries)
        product = np.matmul(rotateMatrix, self.coordinate)
        self.coordinate = product
    
    def rotateAroundX(self, sign):
        entries = [[math.cos(sign * math.pi / 12), 0, math.sin(sign * math.pi / 12), 0],
                    [0, 1, 0, 0],
                    [-math.sin(sign * math.pi / 12), 0, math.cos(sign * math.pi / 12), 0],
                    [0, 0, 0, 1]]
        rotateMatrix = np.array(entries)
        product = np.matmul(rotateMatrix, self.coordinate)
        self.coordinate = product

    def rotateAroundY(self, sign):
        entries = [[math.cos(sign * math.pi / 12), math.sin(sign * math.pi / 12), 0, 0],
                    [0, 0, 1, 0],
                    [-math.sin(sign * math.pi / 12), math.cos(sign * math.pi / 12), 0, 0],
                    [0, 0, 0, 1]]
        rotateMatrix = np.array(entries)
        product = np.matmul(rotateMatrix, self.coordinate)
        self.coordinate = product

def appStarted(app):
    app.atoms = list()
    #            C, H, O, N, P, S
    elemCount = [0, 0, 0, 0, 0, 0]
    app.coordinates, app.elems = coordinatesInNumpy, elemsInNumpy
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
    app.level = -10
    app.isTimerFired = False
    url_home = 'https://mma.prnewswire.com/media/1424831/CRISPeR_Gene_Editing.jpg?w=1200'
    app.image_home = app.loadImage(url_home)
    app.image_home_scaled = app.scaleImage(app.image_home, 2/3)
    url_design_scissors = "https://www.sciencenewsforstudents.org/wp-content/uploads/2019/11/080819_ti_crisprsicklecell_feat-1028x579-1028x579.jpg"
    image_design_scissors = app.loadImage(url_design_scissors)
    # app.organism_image = 
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
    app.dna_bases = {"adenine", "guanine", "cytosine", "thymine"}
    app.rna_bases = {"adenine", "guanine", "cytosine", "uracil"}
    app.amino_acids = {"ala", "arg", "asn", "asp", "cys", "gln", "glu", "gly", 
                       "his", "ile", "leu", "lys", "met", "phe", "pro", "ser",
                       "thr", "trp", "tyr", "val"}
    app.width, app.height = 1000, 800
    app.wantInput = False
    app.timesteps = np.arange(0, 1000, .3)
    resetApp(app)

def resetApp(app):
    # Homepage
    app.isHomepage = True
    app.isDesignPage = False
    app.isHelpPage = False
    app.isIntroPage = False 
    app.drawVisualPage = False
    app.drawCustomSeqPage = False
  
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
    elif event.key == "Left":
        for atom in app.atoms:
            atom.rotateAroundZ(-1)
            atom.coordinate2D = threeDToTwoD(atom.coordinate)
    elif event.key == "Right":
        for atom in app.atoms:
            atom.rotateAroundZ(1)
            atom.coordinate2D = threeDToTwoD(atom.coordinate)
    elif event.key == "r":
        app.isTimerFired = not app.isTimerFired
    elif not app.isHelpPage and event.key == "h":
        app.isHelpPage = True
    elif app.isHelpPage and event.key == "x":
        app.isHelpPage = False
    elif app.isHomepage:
        app.isHomepage = False
        app.isIntroPage = True

# def timerFired(app):
#     while app.isTimerFired:
#         for atom in app.atoms:
#             atom.rotateAroundZ(-0.5)
#             atom.coordinate2D = threeDToTwoD(atom.coordinate)

def getUserInput(app, prompt):
    if app.wantInput:
        return simpledialog.askstring('getUserInput', prompt)
    # app.timerDelay = 10000
    app.wantInput = False

def mousePressed(app, event):
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
            goToDesignPage(app)
        elif inButton(event, app.buttonBottomRight):
            goToVisualPage(app)
    elif app.isVisualPage:
        if inButton(event, app.buttonTopLeft):
            resetApp(app)

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
    elif app.isHelpPage:
        drawHelpPage(app,canvas)

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
    drawCenter1Button(app, canvas, "RoyalBlue4")
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
    drawBottomRightButton(app, canvas, "pale violet red")
    # canvas.create_image(app.width / 3.3, app.height / 1.7,
    #                     image=ImageTk.PhotoImage(app.image_design))
    buttonLeft = app.buttonBottomLeft
    canvas.create_text((buttonLeft[0] + buttonLeft[2])/2, (buttonLeft[1] + buttonLeft[3])/2,
                       text = "Use Template Sequence", fill = "snow", 
                       font = "Arial 15 bold")
    buttonRight = app.buttonBottomRight
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
        y = 100 * np.cos(2 * np.pi * (t + 15) / 255) + app.height/2
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
        # if y < 500 and (t % (1000//3) < (1000//6)):
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
    for atom in app.atoms:
        # Converts 2D coordinates into TK coordinates
        coordinate = (atom.coordinate2D[0] + app.width / 2, 
                      atom.coordinate2D[1] + app.width / 3)
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
    # Displays chemical composition
    counter = 0
    for element in app.percent:
        canvas.create_text(3.5 * app.width / 5, (1 + 0.3 * counter) * app.height / 10, 
                    text = f"{element}: {app.percent[element]}%", anchor = "w", 
                    fill = "mint cream", font = "Arial 15 bold")
        counter += 1

def drawCustomSeqPage(app, canvas):
        # canvas.create_image(3 * app.width / 4 , app.height / 5, 
    #                     image=ImageTk.PhotoImage(app.organism_image))
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
                        text = "Return", font = "Arial 20",fill="snow")

def drawBackground(app, canvas, color):
    canvas.create_rectangle(0, 0, app.width, app.height, fill = f"{color}")

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

# def twoDToTK(app, coordinate):
#     # X remains unchanged
#     x, oldY = coordinate[0], coordinate[1]
#     newY = -oldY + app.height / 2
#     return [x, newY]

def drawAxes(app, canvas):
    # X axis
    canvas.create_line(app.width / 5, 4 * app.height / 5, app.width / 2, app.height / 2)
    # Y axis
    canvas.create_line(app.width / 2, app.height / 2, 4 * app.width / 5, 4 * app.height / 5)
    # Z axis
    canvas.create_line(app.width / 2, app.height / 5, app.width / 2, app.height / 2)

# Atomic radii of the elements
# C =  77, H = 37, O = 73, N = 70, P = 110, S = 103
def drawCarbon(app, canvas, coordinate):
    r = 77/73 + 0.15 * app.level
    (x, y) = (coordinate[0], coordinate[1])
    canvas.create_oval(x-r, y-r, x+r, y+r, fill = 'slategray1')

def drawHydrogen(app, canvas, coordinate):
    r = 37/73 + 0.15 * app.level
    (x, y) = (coordinate[0], coordinate[1])
    canvas.create_oval(x-r, y-r, x+r, y+r, fill = 'floral white')

def drawOxygen(app, canvas, coordinate):
    r = 73/73 + 0.15 * app.level
    (x, y) = (coordinate[0], coordinate[1])
    canvas.create_oval(x-r, y-r, x+r, y+r, fill = 'IndianRed3')

def drawNitrogen(app, canvas, coordinate):
    r = 70/73 + 0.15 * app.level
    (x, y) = (coordinate[0], coordinate[1])
    canvas.create_oval(x-r, y-r, x+r, y+r, fill = 'RoyalBlue3')

def drawPhosphorus(app, canvas, coordinate):
    r = 110/73 + 0.15 * app.level
    (x, y) = (coordinate[0], coordinate[1])
    canvas.create_oval(x-r, y-r, x+r, y+r, fill = 'gold', outline = "gold")

def drawSulfur(app, canvas, coordinate):
    r = 103/73 + 0.15 * app.level
    (x, y) = (coordinate[0], coordinate[1])
    canvas.create_oval(x-r, y-r, x+r, y+r, fill = 'LightPink1', outline = 'LightPink1')

runApp(width = 1000, height = 800)