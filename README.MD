# 15-112 Fundamentals of Programming and Computer Science Term Project
Claire Chen

Video Demo:
https://youtu.be/CYP1r7d8ziw

Project description:
MyCrispr is an app that generates a DNA sequence based on the protein PDB files. With this app, users can select a protein by providing its PDB code, edit the protein's DNA sequence of interest, visualize the editing of the protein-encoding gene, and compare the structure of proteins translated from the sequence before and after the edit. 

Project logistics:
The user should run the file main.py and install the following modules:
- cmu_112_graphics (which is included in the folder)
- numpy
- functools
- PIL
- Bio.PDB
- pathlib

Shortcuts & Commands:
- Homepage: press any key to begin
- Protein selection: use the mouse to click on the button "Choose your protein"
                     and enter a 4-digit PDB protein code when prompted.
- Protein visualization page: use the mouse to click on the button "Visualize protein."
                              Then, either press arrow keys or click on the arrow-shaped
                              buttons to rotate the protein 3-dimensionally. Press the "+"
                              or "-" keys to zoom in or zoom out.
- DNA editing page: after the page loads (after clicking the button "Edit protein"), 
                    press the left or right keys to view other parts of the sequence
                    not shown on the screen. 
- Help page: at all times, press the "h" key to see more directions on use of the app
