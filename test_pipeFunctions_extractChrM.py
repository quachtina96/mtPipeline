import shlex
import subprocess as sp
import os
import numpy
import pipeFunctions as pf

path = "/gpfs/home/quacht/debug/test_getDepth"
depths = []

os.chdir(path)
for subdir, dirs, files in os.walk(path): #go through all the files in the folder, designated by path
	for file in files: #for each file
		if (file.endswith("exome.bam")):
			print "Extracting the chrM region of " + str(file)
			mtExtract=pf.extractChrM(file, path)	
			pf.index(mtExtract)
			print ""
			print "job complete"


