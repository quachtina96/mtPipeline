#python wrapper for the histDepth script
import shlex
import subprocess as sp
import os
import sys
import numpy
import pipelineFunc
import matplotlib

histCommand = "Rscript histDepth.R"
arg = shlex.split(histCommand)
print arg
histJob = sp.Popen(arg,stdout=sp.PIPE)
stdout, stderr = histJob.communicate()
print " job completed "
histogramName = "hist_chrM_" + ".pdf" #need to edit this line to include the name of the file
print histogramName
os.rename("depth.pdf", histogramName)