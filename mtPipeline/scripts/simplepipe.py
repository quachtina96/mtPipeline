#!/usr/bin/env python

# this is the second part of the pipeline, which calculates and outputs the coverage of the chrM portion
# of the exome before extracting it and remapping it to the rCRS.
import shlex
import subprocess as sp
import os
import sys
import numpy
import getopt


def usage():
	print """Merging exome.part.bam files extracting mtDNA and analyzing coverage.
Version 1 Written by Tina Quach, 2015

Options:
	-m      path to mtPipeline directory (include last "/")
	-i      input folder (e.g. ID18_Father)
	-h      view usage
	"""

# ALLOW THE SCRIPT TO TAKE IN PARAMETERS
try:
	opts, args = getopt.getopt(sys.argv[1:], "m:i:")
except getopt.GetoptError, err:
	print str(err)
	usage()
	sys.exit()

# default path
mtPipeDir = "/gpfs/home/quacht/"
sampleDIr = ""

#read in options
for o, a in opts:
	if o == "-h":
		usage()
		sys.exit()
	elif o == "-m":
		mtPipeDir = a
	elif o == "-i":
		sampleDir = a
	else:
		assert False, "Unhandled option."

scriptsDir = str((str(mtPipeDir)+"scripts"))
os.chdir(scriptsDir)
import pipeFunctions as pf

# path to sample directory
path = str(sampleDir)
os.chdir(path)
print "Current directory: " + os.getcwd()
# get the name of sample directory as the name of the sample
sampleName = path.strip().split("/")[-1]
print "Working with sample: " + sampleName
print ""

# MERGE THE PART.BAM FILES
parts = []
# go through all the files in the folder, designated by path
for subdir, dirs, files in os.walk(path):
	for file in files:
		if (file.find(".bai") == -1):
			if (file.find(str(sampleName)) != -1):
				parts.append(str(file))
mergeArgs = " ".join(parts)
mergedBam = sampleName + "_exome.bam"
command = "samtools merge " + mergedBam + \
	" " + mergeArgs
print command
print ""
print "CURRENTLY MERGING THE PART.BAM FILES..."
print ""
args = shlex.split(command)
job = sp.Popen(args, stdout=sp.PIPE)
stdout, stderr = job.communicate()
print sampleName + " part.bam files merged."
print mergedBam + " created."

# GET DEPTH OF chrM REGION & EXTRACT THE MTDNA
for subdir, dirs, files in os.walk(path):
	for file in files:
		if (file.endswith("exome.bam")):
			print "Currently working with " + str(file)
			print "Indexing file..."
			pf.index(file)
			print "file indexed."
			# CALCULATE THE DEPTH (PRE-REMAP)
			depths = pf.getDepths(file)
			stdout = sys.stdout
			depthAnalysis = open(
				str(sampleName + "depthAnalysis_beforeExtract.txt"), "w+")
			sys.stdout = depthAnalysis
			print "Results:"
			print "Before Remapping:"
			pf.analyzeCoverage(file, depths)
			depthAnalysis.close()
			sys.stdout = stdout

			# EXTRACT THE MTDNA mother_exome.bam
			print "Extracting the chrM region of " + str(file)
			mtExtract = pf.extractChrM(file, path)
			pf.index(mtExtract)
			print ""
# REMAP THE EXTRACTED MTDNA TO THE REFRENCE MITOCHONDRIAL GENOME
for subdir, dirs, files in os.walk(path):
	for file in files:  # for each file
		if (file.endswith("exome_mtExtract.bam")):
			print "Remapping " + str(file)
			print str(file)
			csortedRemappedBam = pf.remap2fa(
				str(file), path, reference="/gpfs/home/quacht/data/chrRCRS.fa")
			print csortedRemappedBam + " created"
# GET AND ANALYZE THE DEPTH OF COVERAGE FOR THE REMAPPED CHRM BAM FILES
for subdir, dirs, files in os.walk(path):
	for file in files:
		if (file.endswith("csort.bam")):
			print "Currently working with " + str(file)
			# CALCULATE THE DEPTH (PRE-REMAP)
			depths = pf.getDepths(str(file), extract=True)
			stdout = sys.stdout
			sys.stdout = open(
				str(sampleName + "depthAnalysis_afterExtract.txt"), "w+")
			print "Results:"
			print "After Remapping:"
			pf.analyzeCoverage(file, depths)
			sys.stdout = stdout
			depthAnalysis.close()
