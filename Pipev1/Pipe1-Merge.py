#!/usr/bin/env python

#The first part of the pipeline merges files that fall under the same sample.
#The 3 possible sample types are Mother, Father, and Proband.
#TO USE: set path equal to the path to directory containing only the bam and bai files of interest
# 		 adjust mergeBam variable to match the title you would like for your merged Bam


import shlex
import subprocess as sp
import os
import sys

path = "/gpfs/home/quacht/mergerg"
sp.call('cd /gpfs/home/quacht/mergerg', shell=True) #move into the directory

samples = ["Mother", "Father", "Proband"]
for sample in samples:
	print sample
	parts=[]
	for subdir, dirs, files in os.walk(path): #go through all the files in the folder, designated by path
		for file in files: #for each file
			if (file.find(".bai") == -1): #if the file is a bam file (not the index)
				if (file.find(str(sample)) != -1): 
					parts.append(str(file))
	mergeArgs = " ".join(parts)
	mergedBam = "ID18_" + sample + "_exome.bam"
	command = "samtools merge " + mergedBam + " " + mergeArgs #create the command string
	print command
	stdout,stderr = sp.Popen(shlex.split(command),stdout=sp.PIPE).communicate()
