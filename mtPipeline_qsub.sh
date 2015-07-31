#!/bin/bash

mtPipeline=/gpfs/home/quacht/
mtPipeScripts="${mtPipeline}scripts/"
pathToSampleDirs=/gpfs/home/quacht/ID18/

cd $pathToSampleDirs
pwd
for sampleDir in $(ls); do

echo "Working with $sampleDir"
pathToSampleDir="${pathToSampleDirs}${sampleDir}"
bash "${mtPipeScripts}mtPipeline.sh" -i "${pathToSampleDir}"
pwd

echo "############ ORGANIZING OUTPUTS ##############"
bash "${mtPipeScripts}cleanUp.sh" -i ${pathToSampleDir}
pwd

echo "############### COPY VCF ###############3#"
cd $pathToSampleDirs
if [ ! -d "${pathToSampleDirs}VCF" ] ; then
	 mkdir "${pathToSampleDirs}VCF"
fi

cp *.vcf "${pathToSampleDirs}VCF"

echo "Job complete."
done
