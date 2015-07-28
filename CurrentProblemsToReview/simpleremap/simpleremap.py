#!/usr/bin/env python

# this is the second part of the pipeline, which calculates and outputs the coverage of the chrM portion
# of the exome before extracting it and remapping it to the rCRS.
import shlex
import subprocess as sp
import os
import sys
import numpy
import pipeFunctions as pf

path = "/gpfs/home/quacht/debug"
depths = []

chdir_command = 'cd ' + path
sp.call(chdir_command, shell=True)  # move into the directory
# REMAP
# go through all the files in the folder, designated by path
for subdir, dirs, files in os.walk(path):
    for file in files:  # for each file
        if (file.endswith("Father_exome_mtExtract.bam")):
            print "Remapping " + str(file)
            remappedSam = pf.remap(file, path)
            remappedBam = pf.samtobam(remappedSam)
            print remappedBam
            print ""
            print "getting header"
            header = pf.getHeader(remappedBam)
            print header
            print ""
            print "viewing bam"
            pf.viewBam(remappedBam)
