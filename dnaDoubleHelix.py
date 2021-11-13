# Source: https://matplotlib.org/stable/gallery/user_interfaces/embedding_in_tk_sgskip.html
from cmu_112_graphics import *
# from seqAnalysis import *
# from userInterfaces import *
import tkinter
import numpy as np
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure
from matplotlib import pyplot as plt

app = tkinter.Tk()
app.wm_title("Tk DNA Visualizer (needs to be integrated with cmu_112_graphics)")
strand1 = Figure(figsize=(5, 1), dpi=100)
strand2 = Figure(figsize=(5, 1), dpi=100)
t = np.arange(0, 8, .01)
print(t)
strand1.add_subplot().plot(2 * np.sin(2 * np.pi * t / 3.4))
strand2.add_subplot().plot(2 * np.cos(2 * np.pi * t / 3.4))
canvas1 = FigureCanvasTkAgg(strand1, master=app)
# canvas1.draw()
canvas2 = FigureCanvasTkAgg(strand2, master=app)
# canvas2.draw()
canvas1.get_tk_widget().pack(side=tkinter.TOP, fill=tkinter.BOTH, expand=1)
canvas2.get_tk_widget().pack(side=tkinter.BOTTOM, fill=tkinter.BOTH, expand=1)
tkinter.mainloop()