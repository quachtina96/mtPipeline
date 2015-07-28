#!/usr/bin/env python
"""Given the path to the directory containing the part.bam files (indexed,
 meaning their respective bai files are also included) for 
a specific sample, this script will merge all the files in that folder"""
import shlex
import subprocess as sp
import os
import numpy
import sys
sys.path.insert(0, '/gpfs/home/quacht')
import pipeFunctions as pf

sampleDir = "/gpfs/home/quacht/debug/ID18_Father"
depths = []

os.chdir(sampleDir)
print "Current directory: " + os.getcwd() 
# get the name of sample directory as the name of the sample
sampleName = sampleDir.strip().split("/")[-1]
print sampleName + " is the sample name"
print ""
parts = []
# go through all the files in the folder, designated by path
for subdir, dirs, files in os.walk(sampleDir):
    for file in files:  # for each file
        # if the file is a bam file (not the index)
        if (file.find(".bai") == -1):
            if (file.find(str(sampleName)) != -1):
                parts.append(str(file))
mergeArgs = " ".join(parts)
mergedBam = sampleName + "_exome.bam"
# create the command string
command = "samtools merge " + mergedBam + \
    " " + mergeArgs  
print command
print ""
print "CURRENTLY MERGING THE PART.BAM FILES..."
print ""
stdout, stderr = sp.Popen(command), stdout=sp.PIPE, shell=True).communicate()
print sampleName + " part.bam files merged."
print mergedBam + " created."

