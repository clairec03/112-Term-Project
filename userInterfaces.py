from cmu_112_graphics import *
from main import *

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

def getUserInput(app, prompt):
    return simpledialog.askstring('getUserInput', prompt)

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