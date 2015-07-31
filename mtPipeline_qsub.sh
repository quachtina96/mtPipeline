#!/bin/bash

#This script is the only script that people should need to interact with.
#The variables below should be changed to meet user's needs

#path to the mtPipeline folder
mtPipeline=/gpfs/home/quacht/ 
#path to the mtPipeline scripts
mtPipeScripts="${mtPipeline}scripts/"
#path to the directory that holds the sample directories within (e.g. ID18 holds ID18_Father, ID18_Mother, and ID18_Proband)
pathToSampleDirs=/gpfs/home/quacht/ID18/

cd $pathToSampleDirs
pwd

for sampleDir in $(ls); do
#run the analysis on single sample
echo "Working with $sampleDir"
pathToSampleDir="${pathToSampleDirs}${sampleDir}"
bash "${mtPipeScripts}mtPipeline.sh" -i "${pathToSampleDir}"
pwd

#clean up the sample directory
echo "############ ORGANIZING OUTPUTS ##############"
bash "${mtPipeScripts}cleanUp.sh" -i "${pathToSampleDir}"
pwd

#copy the vcf for that sample to a VCF folder for the sample set for future VCF analysis
echo "############### COPY VCF ###############3#"
cd $pathToSampleDirs
if [ ! -d "${pathToSampleDirs}VCF" ] ; then
	 mkdir "${pathToSampleDirs}VCF"
fi
cp *.vcf "${pathToSampleDirs}VCF"

echo "Job complete."
done
