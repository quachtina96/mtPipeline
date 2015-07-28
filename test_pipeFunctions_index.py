import shlex
import subprocess as sp
import os
import sys
import numpy
import pipeFunctions as pf

path = "/gpfs/home/quacht/debug"
depths = []

os.chdir(path)
for subdir, dirs, files in os.walk(path): #go through all the files in the folder, designated by path
    for file in files: #for each file
        if (file.endswith("csort.bam")):
            print "Currently working with " + str(file)
            pf.index(file)
