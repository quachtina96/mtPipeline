import shlex
import subprocess as sp
import os
import sys
import numpy
from collections import Counter

def remapToChrRCRS(bam, path, reference="chrRCRS.fa"):
    # old working function
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

	
	
#this function calculates the coverage of the mitochondrial region of a bam file and returns 
#an array of the depths (in order of base)
def getDepths(bam, extract=False):
	depthCommand = "samtools depth -r chrM %s" %(bam) #get depth of coverage at chrM contig of the wg              
	if (extract):
		depthCommand = "samtools depth %s" %(bam)
	arg = shlex.split(depthCommand)
	fileOutput = open("depth.txt", 'w+') #open the file that will hold the output of depthCommand
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
		baseCount +=1
		sumDepth += int(fields[2])
		depths.append(int(fields[2]))
	f.close()
	return depths



#this function writes the output of 
def analyzeCoverage(depths):
	print "Depth of Coverage Statistical Analysis"
	#print "average = " + str(sum(depths)/len(depths))  #adds all the depths and divides by number of bases
	print "mean = " + str(numpy.mean(depths))
	print "std = " + str(numpy.std(depths))
	print "minimum = " + str(numpy.min(depths))
	print "max = " + str(numpy.max(depths))
	print "median = " + str(numpy.median(depths))
	data = Counter(depths)
	print "mode = " + str(data.most_common(1))
	histvalues, bins = numpy.histogram(depths, bins=20, range=(0,100))
	print "Histogram Info"
	i = 0
	for y in histvalues:
		print str(bins[i]) + " to " + str(bins[i+1]) + ":	" + str(y)
		i+=1


def extractChrM(bam):
	newName = str(bam[:-4] + "_mtExtract.bam") #create the label for output bam, the mtDNA extract file
	command = "samtools view -b -h " + str(bam) + " chrM " #extract the mtDNA (should be in the sdout)
	newPath = path + "/" + newName
	with open(newPath, 'w+') as outfileHandle:
		process = sp.Popen(shlex.split(command), stdout = outfileHandle).communicate()
	print(str(file)," mtDNA extract complete")
	return newName

def readCount(bam):
	command = "samtools index " + bam #create the command string
	stdout,stderr = sp.Popen(shlex.split(command),stdout=sp.PIPE).communicate()
	command = "samtools view -c " + bam #create the command string
	stdout,stderr = sp.Popen(shlex.split(command),stdout=sp.PIPE).communicate()
	mtReadCount = stdout
	return int(mtReadCount) #make a tuple containing file and stdout

