import shlex
import subprocess as sp
import os
import sys
import numpy
from collections import Counter


def index(bam):
    """ takes in bam file and creates a corresponding bai file."""
    command = "samtools index " + str(file)  # create the command string
    stdout, stderr = sp.Popen(
        shlex.split(command), stdout=sp.PIPE).communicate()


def qsort(bam):
    """takes in bam input and returns name of bam sorted by query name"""
    # bam must be sorted before conversion to fq
    command = "samtools sort -n " + bam + " " + bam[:-4] + ".qsort"
    print "Calling samtools sort..."
    job = shlex.split(command)
    # execute sort command and send output to stdout
    stdout, stderr = sp.Popen(job, stdout=sp.PIPE).communicate()
    return bam[:-4] + ".qsort.bam"


def bamtofastq(bam):
    """ input: qsorted bam, output: 2 fastq files for paired reads"""
    newFastq1 = bam[:-4] + "_1.fq"  # strings holding names of the fastq files
    newFastq2 = bam[:-4] + "_2.fq"
    command = "bedtools bamtofastq -i " + \
        bam[:-4] + ".qsort.bam" + " -fq " + newFastq1 + " -fq2 " + newFastq2
    # ^ convert from bam to fastq
    print "Calling bedtools bamtofastq..."
    job = shlex.split(command)
    # execute sort command and send output to stdout
    stdout, stderr = sp.Popen(job, stdout=sp.PIPE).communicate()
    print "bamtofastq executed"
    return (newFastq1, newFastq2)


def samtobam(sam, reference="chrRCRS.fa"):
    newBam = str(sam[:-4] + ".bam")
    print "Converting sam to bam"
    command = "samtools view -bT " + reference + " " + str(sam) + " -o " + newBam
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

def remap(bam, path, reference="chrRCRS.fa"):
    qsortedBam = qsort(bam)
    fastq_tuple = bamtofastq(qsortedBam)
    command = "bwa mem "+reference + " " + \
        fastq_tuple[0] + " " + fastq_tuple[1]
    print "Calling bwa mem..."
    remappedSam = str(bam[:-4] + "remap.sam")
    args = shlex.split(command)
    newPath = path + "/" + remappedSam
    outfileHandle = open(newPath, 'w+')
    remapJob = sp.Popen(args, stdout=outfileHandle)
    stdout, stderr = remapJob.communicate()
    print "Remap executed."
    return remappedSam


def csort(bam):
    # bam must be sorted before conversion to fq
    command = "samtools sort " + bam + " " + bam[:-4] + ".csort"
    print "Sorting bam"
    job = shlex.split(command)
    # execute sort command and send output to stdout
    stdout, stderr = sp.Popen(job, stdout=sp.PIPE).communicate()
    return bam[:-4] + ".csort.bam"


def getDepths(bam, extract=False):
    """# this function calculates the coverage of the mitochondrial region of a bam file and returns
    # an array of the depths (in order of base)"""
    # get depth of coverage at chrM contig of the wg
    depthCommand = "samtools depth -r chrM %s" % (bam)
    if (extract):
        depthCommand = "samtools depth %s" % (bam)
    arg = shlex.split(depthCommand)
    # open the file that will hold the output of depthCommand
    fileOutput = open("depth.txt", 'w+')
    depthJob = sp.Popen(arg, stdout=fileOutput)  # create job
    depthInitial = depthJob.communicate()  # excecute job
    fileOutput.close()
    # calculate sum of depths
    # will hold array ['chrM', 'currentBase#', 'depth@currentBase']
    fields = []
    sumDepth = 0
    baseCount = 0
    depths = []
    f = open("depth.txt", 'r+')
    # for each position/base in the chrM contig (line in depth file)
    for line in f:
        # split the line into an array of 3 fields
        fields = line.strip().split('\t')
        baseCount += 1
        sumDepth += int(fields[2])
        depths.append(int(fields[2]))
    f.close()
    return depths



def analyzeCoverage(file, depths):
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
        print str(bins[i]) + " to " + str(bins[i+1]) + ":	" + str(y)
        i += 1

    histCommand = "Rscript histDepth.R"
    arg = shlex.split(histCommand)
    histJob = sp.Popen(arg, stdout=sp.PIPE)
    stdout, stderr = histJob.communicate()
    histogramName = "hist_chrM_" + str(file[:-9]) + ".pdf"
    os.rename("depth.pdf", histogramName)
    print "See " + histogramName + "to view histogram regarding the coverage of the chrM region of " + str(file)


def extractChrM(bam):
    # create the label for output bam, the mtDNA extract file
    newName = str(bam[:-4] + "_mtExtract.bam")
    # extract the mtDNA (should be in the sdout)
    command = "samtools view -b -h " + str(bam) + " chrM "
    newPath = path + "/" + newName
    with open(newPath, 'w+') as outfileHandle:
        process = sp.Popen(
            shlex.split(command), stdout=outfileHandle).communicate()
    print(str(file), " mtDNA extract complete")
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
