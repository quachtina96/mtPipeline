#!/usr/bin/env python

import shlex
import subprocess as sp
import os
import sys
import numpy
import pipelineFunc

depths = []
fileList = []
path = "/gpfs/home/quacht/ID18exome/merge"
sp.call('cd /gpfs/home/quacht/ID18exome/merge', shell=True) #move into the directory
for subdir, dirs, files in os.walk(path): #go through all the files in the folder, designated by path
	for file in files: #for each file
		#filepath=os.path.join(subdir,file) #get the file path
		if (file.endswith("exome_mtExtract.bam")):
			print "Remapping " + str(file)
			remappedBamFile = pipelineFunc.remap2mt(file,path)
			print remappedBamFile + "created"
