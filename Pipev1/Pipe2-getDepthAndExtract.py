#!/usr/bin/env python

#this is the second part of the pipeline, which calculates and outputs the coverage of the chrM portion
#of the exome before extracting it and remapping it to the rCRS.
import shlex
import subprocess as sp
import os
import sys
import numpy
import pipelineFunc
import matplotlib

fileList = []
depths = []
path = "/gpfs/home/quacht/ID18exome/merge/chrRCRS_remap"
sp.call('cd /gpfs/home/quacht/ID18exome/merge/chrRCRS_remap', shell=True) #move into the directory
for subdir, dirs, files in os.walk(path): #go through all the files in the folder, designated by path
	for file in files: #for each file
		filepath=os.path.join(subdir,file) #get the file path
		if (file.endswith("exome.bam")):
			print "Currently working with " + str(file)
			print "Indexing file..."
			command = "samtools index " + str(file) #create the command string
			stdout,stderr = sp.Popen(shlex.split(command),stdout=sp.PIPE).communicate()

			#CALCULATE THE DEPTH (PRE-REMAP)
			print "Calculating depth..."
			depthCommand = "samtools depth -r chrM %s" %(file) #get depth of coverage at chrM contig of the wg              
			arg = shlex.split(depthCommand)
			fileOutput = open("depth.txt", 'w+') #open the file that will hold the output of depthCommand
			depthJob = sp.Popen(arg,stdout=fileOutput).communicate() #excecute job
			fileOutput.close()
			#calculate sum of depths
			fields = [] #will hold array ['chrM', 'currentBase#', 'depth@currentBase']
			sumDepth = 0
			baseCount = 0
			depths=[]
			f = open("depth.txt", 'r+')
			for line in f: #for each position/base in the chrM contig (line in depth file)
				fields = line.strip().split('\t') #split the line into an array of 3 fields
				baseCount +=1
				sumDepth += int(fields[2])
				depths.append(int(fields[2]))
			f.close()
			print "Results:"
			print "Before Remapping:"
			pipelineFunc.analyzeCoverage(depths)
			
			histCommand = "Rscript histDepth.R"
			arg = shlex.split(histCommand)
			histJob = sp.Popen(arg,stdout=sp.PIPE)
			stdout, stderr = histJob.communicate()
			histogramName = "hist_chrM_" + str(file[:-9])+ ".pdf"
			os.rename("depth.pdf", histogramName)
			print "See " + histogramName + "to view histogram regarding the coverage of the chrM region of " + str(file)
			

			#EXTRACT THE MTDNA mother_exome.bam
			print "Extracting the chrM region of " + str(file)
			newName = str(file[:-4] + "_mtExtract.bam") #create the label for output bam, the mtDNA extract file
			#print("The new file name is ",newName) #temp
			command = "samtools view -b -h " + str(file) + " chrM " #extract the mtDNA (should be in the sdout)
			newPath = path + "/" + newName
			#print("Now executing ",command, "and writing to ",newPath)
			with open(newPath, 'w+') as outfileHandle:
				process = sp.Popen(shlex.split(command), stdout = outfileHandle).communicate()
			#stdout,stderr = sp.Popen(shlex.split(command),stdout=sp.PIPE).communicate()
			print "Extract complete"
			
			command = "samtools index " + newName #create the command string
			#print(command)
			stdout,stderr = sp.Popen(shlex.split(command),stdout=sp.PIPE).communicate()
			print ""

	






