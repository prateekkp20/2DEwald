# file for running the tests for the approximation of the error functions
import os
from lammps_logfile import File
from pathlib2 import Path 
import numpy as np
import pandas as pd
from random import randint

def replace_line_in_file(file_path, start_words, new_line):
    with open(file_path, 'r') as file:
        lines = file.readlines()

    modified_lines = []
    for line in lines:
        if line.startswith(start_words):
            modified_lines.append(new_line + '\n')
        else:
            modified_lines.append(line)

    with open(file_path, 'w') as file:
        file.writelines(modified_lines)

functionfile = "src/func.cpp"
outputfile = "terminal.txt" 
printfile = "printfile.csv"

for i in range(30,55,1):
    replace_line_in_file(functionfile,"    if(val>","    if(val>"+str(i/10)+"){")
    os.system(f"make clean && make >> {outputfile}")
    os.chdir("run")
    os.system(f"./coulomb.x >> {printfile}")
    os.chdir("..")

os.system(f"rm {outputfile}")