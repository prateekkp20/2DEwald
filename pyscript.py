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
TopHat_File = "tophat.in"
POSCAR_File = "POSCAR_Files/Ewald2D/NaCl/POSCAR.010"

"""Charges"""
NaCl = [1,-1]
CaCl2 = [2,-1]
ChargeFile(ChargeFilePath, NaCl)

"""Output Files"""
CSV_File = "Data/Ewald2D/NaCl/H/H1.csv"

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
# with open(CSV_File, 'wb') as file:
#         for i in range(1,Total):
#             replace_line_in_file(Input_File,"Posfile = ","Posfile = "+POSCAR_File+str(i).zfill(3))
#             for kvecz in range(4,100,1):
#                 replace_line_in_file(Ewald_File,"kz = ","kz = "+str(kvecz))
#                 result = subprocess.run(["./coulomb.x"], stdout=PIPE, stderr=PIPE)
#                 file.write(result.stdout)  # Write stdout in binary format
#                 file.write(result.stderr)  # Write stderr in binary format

"""Experiment D: Polynomial approximation of the erf"""
# with open(CSV_File, 'wb') as file:
#     for i in range(1,Total):
#         for j in range(0,Repeat):
#             replace_line_in_file(Input_File,"Posfile = ","Posfile = "+POSCAR_File+str(i).zfill(3))
#             result = subprocess.run(["./coulomb.x"], stdout=PIPE, stderr=PIPE)
#             if result.returncode != 0:
#                 # print("Error: Integration Library Functions Crashed")
#                 continue
#             file.write(result.stdout)  # Write stdout in binary format
#             file.write(result.stderr)  # Write stderr in binary format

"""Experiment F: Changing the Vacuum size"""
# with open(CSV_File, 'wb') as file:
#     for l in np.arange(55, 100, 0.2):
#         replace_line_in_file(POSCAR_File,"     0.0000000000000000    0.0000000000000000   ","     0.0000000000000000    0.0000000000000000   "+f"{round(l, 12):.12f}"+"0000")
#         result = subprocess.run(["./coulomb.x"], stdout=PIPE, stderr=PIPE)
#         if result.returncode != 0:
#             # print("Error: Integration Library Functions Crashed")
#             continue
#         file.write(result.stdout)  # Write stdout in binary format
#         file.write(result.stderr)  # Write stderr in binary format

"""Experiment G: for some gamma and L, direct 2d ewald corrections"""
# gammavalues = np.arange(0.21,1.001,0.02)
# Lvalues = [380,352,328,310,294,278,266,256,246,236,228,220,214,208,202,198,192,188,184,180,176,174,170,168,164,162,160,156,154,152,150,148,146,146,144,142,140,138,138,136]

# gammavalues = [1,1.5,2,2.5,3,3.5,4,4.5,5,10,20,35,50]
# Lvalues = [140,110,100,90,85,80,75,75,75,70,70,70,70]

# with open(CSV_File, 'wb') as file:
#     # for i in range(1,Total):
#     #     POSCARName = POSCAR_File+str(i).zfill(3)
#     #     replace_line_in_file(Input_File,"Posfile = ","Posfile = "+POSCARName)

#     for gamma, L in zip(gammavalues, Lvalues):
#         replace_line_in_file(TopHat_File,"gamma ",f"gamma = {gamma:.2f}")
#         replace_line_in_file(POSCAR_File,"     0.0000000000000000    0.0000000000000000   ","     0.0000000000000000    0.0000000000000000   "+f"{round(L+10, 12):.12f}"+"0000")

#         for kvecz in range(4,131,1):
#             replace_line_in_file(Ewald_File,"kz = ","kz = "+str(kvecz))
            
#             # result = subprocess.run(["./coulomb.x"], shell=True)
#             result = subprocess.run(["./coulomb.x"], stdout=PIPE, stderr=PIPE)
#             if result.returncode == 0:
#                 file.write(result.stdout)  # Write stdout in binary format
#                 file.write(result.stderr)  # Write stderr in binary format
#                 # print("Error: Integration Library Functions Crashed")
#                 # continue
#             else:
#                 break

"""Experiment H: for some gamma, SPME 2d ewald corrections"""
gammavalues = [0.2, 0.3, 0.4, 0.5, 0.75, 1, 2.5]
with open(CSV_File, 'wb') as file:
    for gamma in gammavalues:
        replace_line_in_file(TopHat_File,"gamma ",f"gamma = {gamma:.2f}")
        
        for orderxy in range(4,16,2):
            replace_line_in_file(Ewald_File,"nx = ","nx = "+str(orderxy))
            replace_line_in_file(Ewald_File,"ny = ","ny = "+str(orderxy))

            for orderz in range(orderxy,16,2):
                replace_line_in_file(Ewald_File,"nz = ","nz = "+str(orderz))

                for gridxy in [2**n for n in range(4, 8)]:
                    replace_line_in_file(Ewald_File,"gx = ","gx = "+str(gridxy))
                    replace_line_in_file(Ewald_File,"gy = ","gy = "+str(gridxy))
                    
                    for gridz in [2**n for n in range(np.log2(gridxy), 12)]:
                        replace_line_in_file(Ewald_File,"gz = ","gz = "+str(gridz))
                        
                        for kvecz in range(4,gridz,1):
                            replace_line_in_file(Ewald_File,"kz = ","kz = "+str(kvecz))

"""Experiment I: Separate positive and negative charges"""