import shlex
import subprocess as sp
import os
import sys
import numpy
import pipeFunctions as pf

path = "/gpfs/home/quacht/debug/test_getDepth"
depths = []

os.chdir(path)
print os.getcwd()
for subdir, dirs, files in os.walk(path): #go through all the files in the folder, designated by path
	for file in files: #for each file
		#filepath=os.path.join(subdir,file) #get the file path
		if (file.endswith("exome_mtExtract.bam")):
			print "Remapping " + str(file)
			print str(file)
			csortedRemappedBam = pf.remap2fa(str(file),path,reference="/gpfs/home/quacht/data/chrRCRS.fa")
			print csortedRemappedBam + " created"