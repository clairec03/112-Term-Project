# Name: Claire Chen
# Andrew ID: ccz
from cmu_112_graphics import *
from seqAnalysis import *
from userInterfaces import *

def appStarted(app):
    url_home = 'https://mma.prnewswire.com/media/1424831/CRISPeR_Gene_Editing.jpg?w=1200'
    app.image_home = app.loadImage(url_home)
    app.image_home_scaled = app.scaleImage(app.image_home, 2/3)
    url_design_scissors = "https://www.sciencenewsforstudents.org/wp-content/uploads/2019/11/080819_ti_crisprsicklecell_feat-1028x579-1028x579.jpg"
    image_design_scissors = app.loadImage(url_design_scissors)
    app.image_design = app.scaleImage(image_design_scissors, 1/3)
    # app.image_design = image_design.crop((200, 200, 400, 579))
    url_pyrimidines = "https://upload.wikimedia.org/wikipedia/commons/thumb/f/f0/Blausen_0324_DNA_Pyrimidines.png/640px-Blausen_0324_DNA_Pyrimidines.png"
    url_purines = "https://upload.wikimedia.org/wikipedia/commons/thumb/9/9f/Blausen_0323_DNA_Purines.png/460px-Blausen_0323_DNA_Purines.png"
    image_cytosine = app.loadImage(url_pyrimidines)
    app.image_cytosine_scaled = image_cytosine.crop((240, 0, 680, 300))
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
    resetApp(app)

def resetApp(app):
    # Homepage
    app.isHomepage = True
    app.isDesignPage = False
    app.isHelpPage = False
    app.isIntroPage = False 

def timerFired(app):
    pass

def keyPressed(app, event):
    if not app.isHelpPage and event.key == "h":
        app.isHelpPage = True
    elif app.isHelpPage and event.key == "x":
        app.isHelpPage = False
    elif app.isHomepage:
        app.isHomepage = False
        app.isIntroPage = True
 
def getUserInput(app, prompt):
    if app.wantInput:
        return simpledialog.askstring('getUserInput', prompt)
    app.wantInput = False

def mousePressed(app, event):
    print(f"event at {event.x}, {event.y}")
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
            print(app.buttonBottomRight)
            print("event is in button on the design page")
            goToPromptWindow(app)
    # if app.isDesignPage:
    #     pass

def goToDesignPage(app):
    app.isIntropage = False
    app.isDesignPage = True

def goToPromptWindow(app):
    app.wantInput = True

def redrawAll(app, canvas):
    drawBackground(app, canvas)
    if app.isHomepage:
        drawHomepage(app, canvas)
    elif app.isDesignPage:
        drawDesignPage(app, canvas)
    elif app.isIntroPage:
        drawIntroPage(app, canvas)
    canvas.create_text(app.width/2, 5 * app.width /6, text = "press h for help",
                    font = "Arial 15",fill="snow")
    if app.isHelpPage:
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
    canvas.create_image(app.width /2 , app.height / 2, 
                        image=ImageTk.PhotoImage(app.image_home))
    button1 = app.buttonCenter1
    canvas.create_text(app.width / 2, 1.8 * app.height / 3,
                    text = "myCrispr", 
                    font = "Helvetica 70", fill="white")
    canvas.create_text(app.width / 2, 2.2 * app.height / 3,
                    text = "CRISPR cas9 enzyme creator at your fingertips", 
                    font = "Arial 16", fill="CadetBlue1")
    canvas.create_rectangle(button1[0], button1[1], button1[2], button1[3],
                            fill = "RoyalBlue4")
    canvas.create_text((button1[0] + button1[2])/2, (button1[1] + button1[3])/2,
                       text="Skip intro", fill="Slategray2", 
                       font = "Arial 15 bold")
    canvas.create_text((button1[0] + button1[2])/2, 0.92 * (button1[1] + button1[3])/2,
                       text="press any key to begin", fill="white",
                       font="Arial 20 bold", anchor = "n")

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
    canvas.create_image(app.width / 3.3, app.height / 1.7,
                        image=ImageTk.PhotoImage(app.image_design))
    buttonLeft = app.buttonBottomLeft
    canvas.create_text((buttonLeft[0] + buttonLeft[2])/2, (buttonLeft[1] + buttonLeft[3])/2,
                       text = "Use Template Sequence", fill = "snow", 
                       font = "Arial 15 bold")
    buttonRight = app.buttonBottomRight
    canvas.create_text((buttonRight[0] + buttonRight[2])/2, (buttonRight[1]+
                        buttonRight[3])/2, text = "Write Custom Sequence", 
                        fill = "snow", font = "Arial 15 bold")

def drawDesignPage(app, canvas):
    drawReturnButton(app,canvas)
    canvas.create_image(app.width / 4, app.height / 4, 
                        image=ImageTk.PhotoImage(app.image_cytosine_scaled))

def drawBottomLeftButton(app, canvas, color):
    button = app.buttonBottomLeft
    canvas.create_rectangle(button[0], button[1], button[2], button[3],
                            fill=f"{color}")

def drawBottomRightButton(app, canvas, color):
    button = app.buttonBottomRight
    canvas.create_rectangle(button[0], button[1], button[2], button[3],
                            fill=f"{color}")

def drawReturnButton(app,canvas):
    buttonReturn = app.buttonTopLeft
    canvas.create_rectangle(buttonReturn[0], buttonReturn[1], buttonReturn[2], 
                            buttonReturn[3], fill = "RoyalBlue3")
    canvas.create_text((buttonReturn[0]+buttonReturn[2])/2, 
                        (buttonReturn[1]+buttonReturn[3])/2,
                        text = "Return", font = "Arial 20",fill="snow")

def drawBackground(app, canvas):
    canvas.create_rectangle(0, 0, app.width, app.height, fill = "mintcream")

runApp(width = 1000, height = 800)