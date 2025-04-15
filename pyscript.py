# file for running the tests
# python3 pyscript.py
import os
from lammps_logfile import File
from pathlib2 import Path 
import numpy as np
import pandas as pd
from random import randint
import sys
import csv
import subprocess
from subprocess import PIPE

def replace_line_in_file(file_path, start_words, new_line):
    """
    Replaces a line in the given file that starts with a specified phrase.
    
    Args:
        file_path (str): Path to the file to be modified.
        start_words (str): The initial words of the line to be replaced.
        new_line (str): The new line to replace the matched line.
    """
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

def ChargeFile(file_path, new_charges):
    """
    Replaces the charges in the given file with new charges.
    
    Args:
        file_path (str): Path to the file to be modified.
        new_charges (list): List of new charges to replace the existing ones.
    """
    with open(file_path, 'w') as file:
        for charge in new_charges:
            file.write(f"{charge}\n")

# Run the make command and check for errors
result = subprocess.run("clear && make clean && make", shell=True)

# Check the return code
if result.returncode != 0:
    print("Error: Makefile compilation failed.")
    sys.exit(1)  # Exit the script with an error code

os.chdir("run/")

"""Input Files"""
ChargeFilePath = "charge.in" 
Input_File  = "input.in"
Ewald_File = "ewald.in"
POSCAR_File = "POSCAR_Files/Ewald2D/NaCl/POSCAR."

"""Charges"""
NaCl = [1,-1]
CaCl2 = [2,-1]
ChargeFile(ChargeFilePath, NaCl)

"""Output Files"""
CSV_File = "Data/Ewald2D/NaCl/D/D1.csv"

"""Constants"""
Total = 80
Repeat = 5

# for i in range(30,55,1):
#     replace_line_in_file(functionfile,"    if(val>","    if(val>"+str(i/10)+"){")
#     os.system(f"make clean && make >> {outputfile}")
#     os.chdir("run")
#     os.system(f"./coulomb.x >> {printfile}")
#     os.chdir("..")

"""Generating the terms for screening function """
# os.system(f"make clean && make >> {outputfile}")
# os.chdir("run")
# os.system(f"./coulomb.x >> {printfile}")
# os.chdir("..")
# os.system(f"rm {outputfile}")

""" run/Data/Ewald2D/NaCl/A/Readme.md """
# with open(CSV_File, 'wb') as file:
#     for i in range(1,Total):
#         replace_line_in_file(Input_File,"Posfile = ","Posfile = "+POSCAR_File+str(i).zfill(3))
#         # for timer in range(0,Repeat):
#         result = subprocess.run(["./coulomb.x"], stdout=PIPE, stderr=PIPE)
#         file.write(result.stdout)  # Write stdout in binary format
#         file.write(result.stderr)  # Write stderr in binary format

"""Checking the convergence of our work with pppm method"""
# with open(CSV_File, 'wb') as file:
#     for orderx in range(8,12,2):
#         replace_line_in_file(Ewald_File,"nx = ","nx = "+str(orderx))

#         # for ordery in range(4,12,2):
#         replace_line_in_file(Ewald_File,"ny = ","ny = "+str(orderx))

#         for orderz in range(4,12,2):
#             replace_line_in_file(Ewald_File,"nz = ","nz = "+str(orderz))

#             for gridx in [2**n for n in range(4, 8)]:
#                 replace_line_in_file(Ewald_File,"gx = ","gx = "+str(gridx))

#                 # for gridy in [2**n for n in range(4, 8)]:
#                 replace_line_in_file(Ewald_File,"gy = ","gy = "+str(gridx))
                
#                 for gridz in [2**n for n in range(4, 8)]:
#                     replace_line_in_file(Ewald_File,"gz = ","gz = "+str(gridz))

                    # for kvecx in range(4,8,1):
                        # replace_line_in_file(Ewald_File,"kx = ","kx = "+str(kvecx))

                        # for kvecy in range(4,8,1):
                            # replace_line_in_file(Ewald_File,"ky = ","ky = "+str(kvecy))

                    # for kvecz in range(4,gridz,1):
                    #     replace_line_in_file(Ewald_File,"kz = ","kz = "+str(kvecz))

                    #     result = subprocess.run(["./coulomb.x"], stdout=PIPE, stderr=PIPE)
                    #     file.write(result.stdout)  # Write stdout in binary format
                    #     file.write(result.stderr)  # Write stderr in binary format

"""Experiment C: Calculating the reciprocal energies using the direct ewald method"""
with open(CSV_File, 'wb') as file:
        for i in range(1,Total):
            replace_line_in_file(Input_File,"Posfile = ","Posfile = "+POSCAR_File+str(i).zfill(3))
            for kvecz in range(4,100,1):
                replace_line_in_file(Ewald_File,"kz = ","kz = "+str(kvecz))
                result = subprocess.run(["./coulomb.x"], stdout=PIPE, stderr=PIPE)
                file.write(result.stdout)  # Write stdout in binary format
                file.write(result.stderr)  # Write stderr in binary format

"""Experiment D: Varying the vacuum in Z direction and"""
with open(CSV_File, 'wb') as file:
    result = subprocess.run(["./coulomb.x"], stdout=PIPE, stderr=PIPE)
    file.write(result.stdout)  # Write stdout in binary format
    file.write(result.stderr)  # Write stderr in binary format