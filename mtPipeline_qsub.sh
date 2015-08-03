#!/bin/bash

#This script is the only script that people should need to interact with.
#The variables below should be changed to meet user's needs
echo pwd
date
#path to the mtPipeline folder
mtPipeline=/gpfs/home/quacht/ 
#path to the mtPipeline scripts
mtPipeScripts="${mtPipeline}scripts/"
#path to the directory that holds the sample directories within (e.g. ID18 holds ID18_Father, ID18_Mother, and ID18_Proband)
pathToSampleDirs=/gpfs/home/quacht/ID18/

bash "${mtPipeScripts}mtPipeline.sh" -i "${pathToSampleDirs}" -m "${mtPipeline}"

date