#!/usr/bin/env python

import shlex
import subprocess as sp
import os
import sys
import numpy
import pipelineFunc


depths = []
fileList = []
path = "/gpfs/home/quacht/ID18exome/merge/chrRCRS_remap"
sp.call('cd /gpfs/home/quacht/ID18exome/merge/chrRCRS_remap', shell=True) #move into the directory
for subdir, dirs, files in os.walk(path): #go through all the files in the folder, designated by path
	for file in files: #for each file
		#filepath=os.path.join(subdir,file) #get the file path
		if (file.endswith("exome_mtExtractremap.csort.bam")):
			#remappedBamFile = pipelineFunc.remap2mt(file,path)
			depthCommand = "samtools depth " + str(file) #get depth of coverage at chrM contig of the wg              
			arg = shlex.split(depthCommand)
			with open("depth.txt", 'w+') as fileOutput: #open the file that will hold the output of depthCommand
				depthJob = sp.Popen(arg,stdout=fileOutput) #create job
			depthInitial = depthJob.communicate() #excecute job
			fileOutput.close()
			#calculate sum of depths
			fields = [] #will hold array ['chrM', 'currentBase#', 'depth@currentBase']
			sumDepth = 0
			baseCount = 0
			depths=[]
			f = open("depth.txt", 'r+')
			for line in f: #for each position/base in the chrM contig (line in depth file)
				fields = line.strip().split('\t') #split the line into an array of 3 fields
				depths.append(int(fields[2]))
			f.close()
			fileList.append(file)
			
			print "Results:"
			print "After Remapping: "
			pipelineFunc.analyzeCoverage(depths)
			
			histCommand = "Rscript histDepth.R"
			arg = shlex.split(histCommand)
			histJob = sp.Popen(arg,stdout=sp.PIPE)
			stdout, stderr = histJob.communicate()
			histogramName = "hist_chrM_" + str(file[:-12])+ ".pdf"
			os.rename("depth.pdf", histogramName)
			print "See " + histogramName + "to view histogram regarding the coverage for " + str(file)






