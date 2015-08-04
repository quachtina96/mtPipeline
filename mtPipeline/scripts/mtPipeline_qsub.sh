#!/bin/bash

#This script is the only script that people should need to interact with.
#The variables below should be changed to meet user's needs
echo pwd
date

#path to parameters
param=/gpfs/home/quacht/scripts/parameters.sh
source $param
#path to the directory that holds the sample directories within (e.g. ID18 holds ID18_Father, ID18_Mother, and ID18_Proband)
pathToSampleDirs=/gpfs/home/quacht/ID18/
#run mtPipeline
bash "${mtPipelineScripts}mtPipeline.sh" -i "${pathToSampleDirs}" -p "${param}" >> log.txt

date