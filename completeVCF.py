import sys
import os
import glob
import math
import re
import ast
from collections import OrderedDict
import getopt

def usage():
    print """This script reads a combined.vcf file (which displays "." or "./." to indicate missing information) and 
    creates a completed.vcf file that includes all information on each sample at positions of interest.  
	It calculates this info from info extracted from the mtDNAassembly-table for the corresponding sample. 

	The methods used to perform calculations were adapted from those in mtVariantCaller.py.

	The script has been tested only with sample trios (father, mother, proband).

	 Written by Tina Quach - 2015		
		
		Mandatory Input Options:
		-c path to combined VCF 
		-t path to dir containing only each sample's mtDNAassembly-table.txt
		"""

#take in arguments
try:
    opts, args = getopt.getopt(sys.argv[1:], "h:c:t:p:")
except getopt.GetoptError, err:
    print str(err)
    usage()
    sys.exit()

combinedVCF = ""
assembleDir = ""

for o, a in opts:
    if o == "-c":
        combinedVCF = a
    if o == "-t":
        assembleDir = a
if (assembleDir == ""):
    print "Missing Required Option: -t <path to dir containing only each sample's mtDNAassembly-table.txt>"
    sys.exit()
if (combinedVCF == ""):
    print "Missing Required Option: -c <path to combined VCF>"
    sys.exit()

# FUNCTIONS ####################################################

# Wilson confidence interval lower bound
def CIW_LOW(het, covBase):
    '''The function calculates the heteroplasmic fraction and the related
    confidence interval with 95% of coverage probability,
    considering a Wilson score interval when n<=40 
    CIw= [1/(1+(1/n)*z^2)] * [p + (1/2n)*z^2 +- z(1/n *(p*q) + ((1/(4n^2))*z^2))^1/2]
    '''
    p = float(het)
    n = float(covBase)
    z = 1.96
    q = 1-het
    num = p*q
    squarez = z*z
    squaren = n*n
    wilsonci_low = round(
        (p+(z*z)/(2*n)-z*(math.sqrt(p*q/n+(z*z)/(4*(n*n)))))/(1+z*z/n), 3)
    if wilsonci_low < 0.0:
        return 0.0
    else:
        return wilsonci_low

# Wilson confidence interval upper bound
def CIW_UP(het, covBase):
    '''The function calculates the heteroplasmic fraction and the related
    confidence interval with 95% of coverage probability,
    considering a Wilson score interval when n<=40 
    CIw= [1/(1+(1/n)*z^2)] * [p + (1/2n)*z^2 +- z(1/n *(p*q) + ((1/(4n^2))*z^2))^1/2]
    '''
    p = het
    n = covBase
    z = 1.96
    q = 1-het
    num = p*q
    squarez = z*z
    squaren = n*n
    wilsonci_up = round(
        (p+(z*z)/(2*n)+z*(math.sqrt(p*q/n+(z*z)/(4*(n*n)))))/(1+z*z/n), 3)
    if wilsonci_up > 1.0:
        return 1.0
    else:
        return wilsonci_up

# Agresti-Coull confidence interval lower bound
def CIAC_LOW(cov, covBase):
    '''The function calculates the heteroplasmic fraction and the related confidence interval 
    for heteroplasmic fraction with 95% of coverage probability,considering the 
    Agresti-Coull interval when n>40'''
    z = 1.96
    n = int(covBase)
    X = float(cov)+(z*z)/2
    # print X, "X"
    N = n+(z*z)
    # print N, "N"
    P = X/N
    # print P, "P"
    Q = 1-P
    # print Q, "Q"
    agresticoull_low = round(P-(z*(math.sqrt(P*Q/N))), 3)
    if agresticoull_low < 0.0:
        return 0.0
    else:
        return agresticoull_low

# Agresti-Coull confidence interval upper bound
def CIAC_UP(cov, covBase):
    '''The function calculates the heteroplasmic fraction and the related confidence interval 
    for heteroplasmic fraction with 95% of coverage probability,considering the 
    Agresti-Coull interval when n>40'''
    z = 1.96
    n = float(covBase)
    X = float(cov)+(z*z)/float(2)
    # print "X",X
    N = n+(z*z)
    # print "N", N
    P = X/N
    # print "P", P
    Q = 1-P
    # print "Q",Q
    agresticoull_up = round(P+(z*(math.sqrt(P*Q/N))), 3)
    if agresticoull_up > 1.0:
        return 1.0
    else:
        return agresticoull_up


# Heteroplasmic fraction quantification
def heteroplasmy(cov, covBase):
    try:
        if covBase >= cov:
            Heteroplasmy = float(cov)/float(covBase)
            het = round(Heteroplasmy, 3)
            return het
        else:
            return 1.0
    except ZeroDivisionError:
        het = 1.0
        return het


def getCovVar(varBase, baseCounts):
    variantDict = {}
    variantDict["A"] = baseCounts[0]
    variantDict["C"] = baseCounts[1]
    variantDict["T"] = baseCounts[2]
    variantDict["G"] = baseCounts[3]

    return variantDict[str(varBase)]


def getVariantInfo(sampleTable, pos):
    """This function reads info in from the assembly table and '
    calculates the necessary info to augment the VCF

    INPUT: position of interest (e.g. 222) and sample assemblyTable 
    for particular sample

    OUTPUT: dictionary whose keys are 
    GT = Genotype (value is int)
    DP = Depth of Coverage (value is int)
    HF = Heteroplasmy Frequency (float)
    CIUP = Upper limit of confidence interval (float)
    CILOW = Lower limit of confidence interval (float)"""

    posInfo = sampleTable[int(pos)]
    # print posInfo
    varBase = str(posInfo).split("\t")[2]  # get the alternate base
    # get depth of coverage at that position
    covBase = int(str(posInfo).split("\t")[3])
    baseCounts = str(posInfo).split("\t")[5].strip()
    baseCounts = baseCounts[1:-1]  # removing the parentheses(38, 0, 0, 0)
    baseCounts = re.sub(",", "", baseCounts)
    baseCounts = baseCounts.split()
    # number of reads supporting variant
    covVar = getCovVar(varBase, baseCounts)
    covTuple = (covVar, covBase)
    hetfreq = heteroplasmy(covVar, covBase)
    if covBase <= 40:
        het_ci_low = CIW_LOW(hetfreq, covBase)
        het_ci_up = CIW_UP(hetfreq, covBase)
    else:
        het_ci_low = CIAC_LOW(covVar, covBase)
        het_ci_up = CIAC_UP(covVar, covBase)

    VarInfo = {}
    VarInfo["GT"] = 0
    VarInfo["DP"] = covBase
    VarInfo["HF"] = hetfreq
    VarInfo["CIUP"] = het_ci_up
    VarInfo["CILOW"] = het_ci_low
    return VarInfo


def addVariantInfo(VarInfo, FieldOrder, VCFcolumn, row):
    """INPUTS: variant info dictionary generated by getVariantInfo,
                    string that tells order of fields desired (e.g."GT:CILOW:CIUP:DP:HF"), 
                    int representing column in VCF that will recieve new variant info,
                    an array holding each element of the line in VCF"""
    infoString = ""
    for x in FieldOrder.split(":"):
        if(infoString == ""):
            infoString = str(VarInfo[x])
        else:
            infoString = infoString + ":" + str(VarInfo[x])
    print "infoString: %s" % (infoString)
    row[VCFcolumn] = infoString
    return row

def readVCF(pathtoVCF):
    """takes in a path to a VCF and returns an array in which each element is a line in the VCF file.
    """
    pathtoVCFArray = pathtoVCF.split("/")
    VCFfileName = pathtoVCFArray.pop()
    print VCFfileName
    VCFdir = "/".join(pathtoVCFArray) + "/"
    print VCFdir
    savedDir = os.getcwd()
    os.chdir(VCFdir)
    VCFfile = open(VCFfileName, 'r')
    VCFlines = []
    for line in VCFfile:
        VCFlines.append(line)
    VCFfile.close()
    os.chdir(savedDir)
    
    return VCFlines

def switchSampleOrder(pathToVCF, currentSampleOrder, desiredSampleOrder):
    """This function takes in 
     the path to a VCF file
     a string containing the current sample order (comma separated)
     and a string containing the desired sample order (comma separated)
    """
    changedVCF=[]
    currList=currentSampleOrder.split(",")
    print currList
    desiredList=desiredSampleOrder.split(",")
    print desiredList
    i=0
    while (i < len(currList)):
        print "ON PAIR %s" %(i+1)
        print "current list is " + currList[i] + " at " + str(i)
        print "desired list is " + desiredList[i] + " at " + str(i) 
        if (currList[i]==desiredList[i]):
            print ""
            print "currlist[i] is equal to desiredList[i]; VCF at sample column %s matches already." %(i+1)     
        else:       
            VCFlines=readVCF(pathToVCF)
            for line in VCFlines:
                if (line.startswith("##")):
                    changedVCF.append(line) 
                    print "header line; skip"               
                else:
                    print "currently on line"
                    print line
                    row = line.strip().split("\t")
                    print "changing sample %s to %s" %(i+1, desiredList[i])
                    currIndex=0
                    while currIndex < len(currList):
                        if currList[currIndex] == desiredList[i]:
                            temp=row[9+i]   
                            row[9+i]=row[currIndex+9]
                            row[currIndex+9]=temp
                            newLine="\t".join(row)
                            changedVCF.append(newLine)  
                            print "NEW LINE:"
                            print newLine
                            temp=currList[i]
                            currList[currIndex]=temp
                            currList[i]=currList[currIndex]
                        currIndex=currIndex+1               
        i=i+1
    newVCF=open("sampleOrderChanged.vcf", "w+")
    for line in changedVCF:
        if (line.endswith("\n")):
            newVCF.write(line)
        else:
            newVCF.write(line +"\n")
    newVCF.close()

#################################################### SCRIPT ##############

# READ IN FILES
#read in tables-c
os.chdir(assembleDir)
tableArray = os.listdir(assembleDir)

motherTable = []
fatherTable = []
probandTable = []
for table in tableArray:
    print table
    f = open(table, 'r')
    if (table.find("Mother") != -1):
        for line in f:
            motherTable.append(line)
        print "motherTable created "
    if (table.find("Father") != -1):
        for line in f:
            fatherTable.append(line)
        print "fatherTable created "
    if (table.find("Proband") != -1):
        for line in f:
            probandTable.append(line)
        print "probandTable created "
    f.close()
print ""

#read in VCF
savedDir = os.getcwd()
pathtoVCFArray = combinedVCF.split("/")
VCFfileName = pathtoVCFArray.pop()
print VCFfileName
VCFdir = "/".join(pathtoVCFArray) + "/"
print VCFdir


os.chdir(VCFdir)
VCFfile = open(VCFfileName, 'r')
VCFlines = []
for line in VCFfile:

    VCFlines.append(line)
VCFfile.close()
print "Printing VCFlines"

# go through each line in the VCF and create a new version that will be
# written to a VCF file
positions = []
completeVCF = []

for line in VCFlines:
    father = False
    mother = False
    proband = False
    if (line.startswith("#")):
        completeVCF.append(line)
    else:
        print "THIS IS THE CURRENT LINE "
        print line
        # create an array of elements from the line
        row = line.strip().split("\t")
        pos = row[1]  # get the start position for variant
        positions.append(pos)
        if (row[9] == "." or row[9] == "./."):  # father
            print "Adding variant info for father for position " + pos
            varInfoDict = getVariantInfo(fatherTable, pos)
            newRow = addVariantInfo(varInfoDict, "GT:DP:CILOW:CIUP:HF", 9, row)
            father = True
        if (row[10] == "." or row[10] == "./."):  # mother
            print "Adding variant info for mother for position " + pos
            varInfoDict = getVariantInfo(motherTable, pos)
            # if row has been updated already
            if(father):
                print "father has been updated already so using newRow"
            # update the new row, not the old one
                newRow = addVariantInfo(
                    varInfoDict, "GT:DP:CILOW:CIUP:HF", 10, newRow)
                print "new row: "
                print newRow
            else:
                print "first sample to be fixed"
                newRow = addVariantInfo(
                    varInfoDict, "GT:DP:CILOW:CIUP:HF", 10, row)
                print "new row: "
                print newRow
            mother = True
        if (row[11] == "." or row[11] == "./."):  # proband
            print "Adding variant info for proband for position " + pos
            varInfoDict = getVariantInfo(probandTable, pos)
            if(father or mother):
                print "using new"
                newRow = addVariantInfo(
                    varInfoDict, "GT:DP:CILOW:CIUP:HF", 11, newRow)
                print "new row: "
                print newRow
            else:
                print "first sample to be fixed"
                newRow = addVariantInfo(
                    varInfoDict, "GT:DP:CILOW:CIUP:HF", 11, row)
                print "new row: "
                print newRow
            proband = True
        if (proband or mother or father):
            newLine = "\t".join(newRow)
            print ""
            print "appending new line to completeVCF"
            print newLine
            completeVCF.append(newLine)
            print ""
            print ""
            print ""
        else:
        	print "appending old line to completeVCF"
        	completeVCF.append(line)
#write complete VCF into a new file
completeVCFfile=open("completed.vcf", "w+")
for line in completeVCF:
	if line.startswith("#"):
		completeVCFfile.write(line)
	elif (line.endswith("\n")):
		completeVCFfile.write(line)
	else:
		completeVCFfile.write(line +"\n")
completeVCFfile.close()

print "completed.vcf created."
print os.getcwd()

orderedVCF=[]
for line in completeVCF:
    if line.startswith("##"):
        orderedVCF.append(line)
    else:
        row = line.strip().split("\t")
        temp=row[9]   
        row[9]=row[11]
        row[11]=temp
        newLine="\t".join(row)
        orderedVCF.append(newLine) 
orderedVCFfile= open(VCFdir+"ordered.vcf", "w+")
for line in orderedVCF:
    if (line.endswith("\n")):
        orderedVCFfile.write(line)
    else:
        orderedVCFfile.write(line +"\n")
orderedVCFfile.close()
os.chdir(savedDir)