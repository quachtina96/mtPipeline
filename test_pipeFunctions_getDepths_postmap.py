import shlex
import subprocess as sp
import os
import numpy
import sys
import pipeFunctions as pf

path = "/gpfs/home/quacht/debug/test_getDepth"
depths = []

os.chdir(path)
for subdir, dirs, files in os.walk(path): #go through all the files in the folder, designated by path
	for file in files: #for each file
		if (file.endswith("csort.bam")):
			print "Currently working with " + str(file)
			#CALCULATE THE DEPTH (PRE-REMAP)
			depths = pf.getDepths(str(file), extract=True)
			print "Results:"
			print "After Remapping:"
			pf.analyzeCoverage(file, depths)


