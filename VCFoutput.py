#!/usr/bin/env python

import getopt, os, ast, sys
from mtVariantCaller import VCFoutput

def usage():
	print """Produces VCF file output from VCF_dict_tmp file. 
		Version 1.1 - Written by Domenico Simone and Claudia Calabrese - 2013-2014
		Version 1.2 - Written by Tina Quach - 2015, renames the VCF based on the sample name		
		Options:
		-s      sample name
		"""

reference_sequence="RCRS"
sampleName=""

try:
	opts, args = getopt.getopt(sys.argv[1:], "h:r:s:")
except getopt.GetoptError, err:
	print str(err)
	usage()
	sys.exit()

for o,a in opts:
	if o =="-s":
		sampleName = a
	else:
		print "Sample Name is required"
		sys.exit()


VCF_dict = ast.literal_eval(open('VCF_dict_tmp', 'r').read())
VCFoutput(VCF_dict, reference=reference_sequence)
currDir= os.getcwd()
oldVCF = "\VCF_file.vcf"
newVCF = sampleName+"\_VCF.vcf"
os.rename(currDir+oldVCF, currDir+newVCF)
