#!/bin/bash

mtPipeline=/gpfs/home/quacht/
mtPipeScripts="${mtPipeline}scripts/"
pathToSampleDirs=/gpfs/home/quacht/ID18/


for sampleDir in $(ls); do

pathToSampleDir="${pathToSampleDirs}${sampleDir}"
bash "${mtPipeScripts}mtPipeline.sh" -i "${pathToSampleDirs}${sampleDir}"

echo "############ ORGANIZING OUTPUTS ##############"

cd ${pathToSampleDir}
partBamDir="${pathToSampleDir}/partbam/"
mkdir $partBamDir
mv *part.bam $partBamDir
mv *part.bam.bai $partBamDir
echo "/partbam holds part.bam files and their indexes"

echo ""
mtExtract="${pathToSampleDir}/mtExtract/"
mkDir $mtExtract
mv *mtExtract* $mtExtract
mv *mtExtractremap.pdf ${pathToSampleDir}
echo "/mtExtract holds the bam/sam/fastq files resulting from coverage analysis"

echo ""
myMtoolbox_out="$pathToSampleDir/myMtoolbox_out"
mkdir $myMtoolbox_out
mv *.rg.* $myMtoolbox_out 
echo "/myMtoolbox_out holds the bam/sam outputs of the indel realignment and marking of duplicates"

echo ""









done
