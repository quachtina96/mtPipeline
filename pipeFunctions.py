import shlex
import subprocess as sp
import os
import sys
import numpy
from collections import Counter
import shutil


def index(bam):
	""" takes in string containing name of bam file (sorted by coordinate) 
	and creates a corresponding bai file."""
	command = "samtools index " + str(bam)  # create the command string
	print command
	args= shlex.split(command)
	job =sp.Popen(args, stdout=sp.PIPE)
	stdout, stderr = job.communicate()

def qsort(bam):
    command = "samtools sort -n " + str(bam) + " " + str(bam)[:-4] + ".qsort"
    print "Calling samtools sort..."
    job = shlex.split(command)
    # execute sort command and send output to stdout
    stdout, stderr = sp.Popen(job, stdout=sp.PIPE).communicate()
    qsortedBam = bam[:-4] + ".qsort.bam"
    return qsortedBam


def bamtofastq(bam):
	""" input: qsorted bam, output: 2 fastq files for paired reads"""
	newFastq1 = str(bam)[:-4] + "_1.fq"  # strings holding names of the fastq files
	newFastq2 = str(bam)[:-4] + "_2.fq"
	command = "bedtools bamtofastq -i " + \
		str(bam)[:-4] + ".qsort.bam" + " -fq " + newFastq1 + " -fq2 " + newFastq2
	# ^ convert from bam to fastq
	print "Calling bedtools bamtofastq..."
	job = shlex.split(command)
	# execute sort command and send output to stdout
	stdout, stderr = sp.Popen(job, stdout=sp.PIPE).communicate()
	print "bamtofastq executed"
	return (newFastq1, newFastq2)


def samtobam(sam, reference="chrRCRS.fa"):
	"""doesn't seem to be working properly"""
	newBam = str(sam[:-4] + ".bam")
	print "Converting sam to bam"
	command = "samtools view -bT " + reference + \
		" " + str(sam) + " -o " + newBam
	args = shlex.split(command)
	job = sp.Popen(args, stdout=sp.PIPE, stderr=sp.PIPE)
	stdout, stderr = job.communicate()
	return newBam


def getHeader(bam):
	command = "samtools view -H " + bam
	args = shlex.split(command)
	stdout, stderr = sp.Popen(args, stdout=sp.PIPE).communicate()
	print "This is the header of " + str(bam)
	print stdout


def viewBam(bam):
	command = "samtools view " + bam
	args = shlex.split(command)
	stdout, stderr = sp.Popen(args, stdout=sp.PIPE).communicate()
	print stdout


def csort(bam):
	# bam must be sorted before conversion to fq
	command = "samtools sort " + bam + " " + bam[:-4] + ".csort"
	print "Sorting bam"
	job = shlex.split(command)
	# execute sort command and send output to stdout
	stdout, stderr = sp.Popen(job, stdout=sp.PIPE).communicate()
	return bam[:-4] + ".csort.bam"


def remapToFasta(bam, path, pathtoreference="/gpfs/home/quacht/data/chrRCRS.fa"):
	"""This function is not working. Please ignore."""
	qsortedBam = qsort(bam)
	fastq_tuple = bamtofastq(qsortedBam)
	command = "bwa mem " + pathtoreference + " " + fastq_tuple[0] + " " + fastq_tuple[1]
	print "Calling bwa mem..."
	remappedSam = str(bam[:-4] + "remap.sam")
	args = shlex.split(command)
	newPath = path + "/" + remappedSam
	outfileHandle = open(newPath, 'w+')
	remapJob = sp.Popen(args, stdout=outfileHandle)
	stdout, stderr = remapJob.communicate()
	print "Remap executed."

	remappedBam = str(bam[:-4] + "remap.bam")
	print "Converting remapped sam to bam"
	command = "samtools view -bT " + pathtoreference + " " + remappedSam + " -o " + remappedBam
	args = shlex.split(command)
	job = sp.Popen(args, stdout=sp.PIPE, stderr=sp.PIPE)
	stdout, stderr = job.communicate()

	csortedBam = csort(remappedBam)
	return csortedBam

def remap2fa(bam, path, reference="/gpfs/home/quacht/data/chrRCRS.fa"):
    """This is the function remaps a given bam file to a reference fasta file. 
    THREE INPUTS: a string containing bam file name, string containing
     path to the sample file folder, and string containing path to the reference fasta.
     OUTPUT: string containing name of bam file remapped to reference and sorted by coordinate."""

    print "Remapping the extracted chrM to rCRS..."
    # remap to the mtGenome (fasta file)
    newFastq1 = bam[:-4] + "_1.fq"  # strings holding names of the fastq files
    newFastq2 = bam[:-4] + "_2.fq"
    # bam must be sorted before conversion to fq
    command = "samtools sort -n " + bam + " " + bam[:-4] + ".qsort"
    print "Calling samtools sort..."
    job = shlex.split(command)
    # execute sort command and send output to stdout
    stdout, stderr = sp.Popen(job, stdout=sp.PIPE).communicate()

    command = "bedtools bamtofastq -i " + \
        bam[:-4] + ".qsort.bam" + " -fq " + newFastq1 + " -fq2 " + newFastq2
    # ^ convert from bam to fastq
    print "Calling bedtools bamtofastq..."
    job = shlex.split(command)
    # execute sort command and send output to stdout
    stdout, stderr = sp.Popen(job, stdout=sp.PIPE).communicate()
    # remap the mtDNA extract to the reference genome
    command = "bwa mem " + reference + " " + newFastq1 + " " + newFastq2
    print "Calling bwa mem..."
    remappedSam = str(bam[:-4] + "remap.sam")
    args = shlex.split(command)
    newPath = path + "/" + remappedSam
    outfileHandle = open(newPath, 'w+')
    remapJob = sp.Popen(args, stdout=outfileHandle)
    stdout, stderr = remapJob.communicate()
    print "Remap executed."

    remappedBam = str(bam[:-4] + "remap.bam")
    print "Converting remapped sam to bam"
    command = "samtools view -bT " + reference + " " + remappedSam + " -o " + remappedBam
    args = shlex.split(command)
    job = sp.Popen(args, stdout=sp.PIPE, stderr=sp.PIPE)
    stdout, stderr = job.communicate()

    # bam must be sorted before conversion to fq
    command = "samtools sort " + remappedBam + \
        " " + remappedBam[:-4] + ".csort"
    print "Sorting remapped bam"
    job = shlex.split(command)
    # execute sort command and send output to stdout
    stdout, stderr = sp.Popen(job, stdout=sp.PIPE).communicate()

    return remappedBam[:-4] + ".csort.bam"
def getDepths(bam, extract=False):
	"""# this function calculates the coverage of the mitochondrial region of a bam file and returns
	# an array of the depths (in order of base)"""

	# get depth of coverage at chrM contig of the wg
	depthCommand = "samtools depth -r chrM %s" % (str(bam))
	if (extract):
		depthCommand = "samtools depth %s" % (str(bam))
	print depthCommand
	arg = shlex.split(depthCommand)
	# open the file that will hold the output of depthCommand
	with open("depth.txt", 'w+') as fileOutput:
		depthJob = sp.Popen(arg, stdout=fileOutput)  # create job
		depthInitial = depthJob.communicate()  # excecute job
	# calculate sum of depths
	# which will hold array ['chrM', 'currentBase#', 'depth@currentBase']
	fields = []
	sumDepth = 0
	baseCount = 0
	depths = []
	f = open("depth.txt", 'r+')
	# for each position/base in the chrM contig (line in depth file)
	for line in f:
		# split the line into an array of 3 fields
		fields = line.strip().split('\t')
		depths.append(int(fields[2]))
	f.close()
	return depths


def analyzeCoverage(file, depths):
	"""This function runs a statistical analysis of a list of depths. 
	INPUT: bam file, list of depths (in order of corresponding base (e.g. depths[0] is the depth at base 1))
	OUTPUT: prints mean, std, minimum, max, median, mode; creates histogram pdf regarding the coverage of the chrM
	region of the inputted bam file"""
	print "Depth of Coverage Statistical Analysis"
	# print "average = " + str(sum(depths)/len(depths))  #adds all the depths
	# and divides by number of bases
	print "mean = " + str(numpy.mean(depths))
	print "std = " + str(numpy.std(depths))
	print "minimum = " + str(numpy.min(depths))
	print "max = " + str(numpy.max(depths))
	print "median = " + str(numpy.median(depths))
	data = Counter(depths)
	print "mode = " + str(data.most_common(1))
	histvalues, bins = numpy.histogram(depths, bins=20, range=(0, 100))
	print "Histogram Info"
	i = 0
	for y in histvalues:
		print str(bins[i]) + " to " + str(bins[i+1]) + ":   " + str(y)
		i += 1

	histCommand = "Rscript /gpfs/home/quacht/scripts/histDepth.R"
	arg = shlex.split(histCommand)
	histJob = sp.Popen(arg, stdout=sp.PIPE)
	stdout, stderr = histJob.communicate()
	histogramName = "hist_chrM_" + str(file[:-10]) + ".pdf"
	os.rename("depth.pdf", histogramName)
	print "See " + histogramName + " to view histogram regarding the coverage of the chrM region of " + str(file)


def extractChrM(bam,path):
	"""This function extracts the chrM region of an exome.bam file
	INPUT: string containing bam file name, path to the directory containing the bam file (should be the sample's directory)
	OUTPUT: the extracted chrM bam file, which ends in '_mtExtract.bam' """
	# create the label for output bam, the mtDNA extract file
	newName = str(bam[:-4] + "_mtExtract.bam")
	# extract the mtDNA (should be in the sdout)
	command = "samtools view -b -h " + str(bam) + " chrM "
	newPath = path + "/" + newName
	with open(newPath, 'w+') as outfileHandle:
		process = sp.Popen(
			shlex.split(command), stdout=outfileHandle).communicate()
	print(str(bam), " mtDNA extract complete")
	return newName


def readCount(bam):
	command = "samtools index " + bam  # create the command string
	stdout, stderr = sp.Popen(
		shlex.split(command), stdout=sp.PIPE).communicate()
	command = "samtools view -c " + bam  # create the command string
	stdout, stderr = sp.Popen(
		shlex.split(command), stdout=sp.PIPE).communicate()
	mtReadCount = stdout
	return int(mtReadCount)  # make a tuple containing file and stdout

def clean(PathToSampleDirectory):
	"""This function removes the files the side products of remapping 
	to ref mitochondrial genome"""
	if (PathToSampleDirectory[-1] != "/" ):
		sideProducts = PathToSampleDirectory + "/temp"
		results = PathToSampleDirectory + "/results"
	else:
		sideProducts = PathToSampleDirectory + "temp"
		results = PathToSampleDirectory + "results"
	if not (os.path.isdir(sideProducts)):
			os.mkdir(sideProducts)
	filesToRemove = []
	for subdir, dirs, files in os.walk(PathToSampleDirectory): 
		for file in files:
			if (file.endswith("qsort.bam") \
				or file.endswith("mtExtract_1.fq") \
				or file.endswith("mtExtract_2.fq") \
				or file.endswith("remap.bam") \
				or file.endswith("remap.sam")):
				newFilePath = sideProducts + "/" + file
				if not os.path.exists(newFilePath):
					shutil.move(PathToSampleDirectory + "/" + file,sideProducts)
					print file
			if (file.endswith("txt") \
				or file.endswith("pdf") \
				or file.endswith("csort.bam")):
				shutil.move(PathToSampleDirectory + "/" + file, results)	
	shutil.rmtree(sideProducts)
	print "removed from %s" %(PathToSampleDirectory)

def setupTest(PathToSampleDirectory):
	if not (os.getcwd() == PathToSampleDirectory):
		os.chdir(PathToSampleDirectory)
	for subdir, dirs, files in os.walk(PathToSampleDirectory): 
		for file in files:
			if not (file.endswith("part.bam") or \
			 file.endswith("part.bam.bai") or \
			 file.endswith(".py") or\
			 file.endswith(".R") or
			 file.endswith (".sh")):
			 os.remove(file)





				

