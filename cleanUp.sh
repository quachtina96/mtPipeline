#!/bin/bash

set -e
set -o pipefail

usage()
{
USAGE="""This script takes in the path to a sample directory and organizes the 'side products'
		of mtPipeline.

		-i option REQUIRED.

		Four directories are created in your sample directory:
		1) partbam (contains all the part.bam and part.bam.bai files)
		2) mtExtract (contains all the files created (excluding results) during the initial coverage analysis in simplepipe.py)
		3) myMtoolbox_out (contains all the files created through myMtoolbox.sh except for results)
		4) results (contains all relevant results to be seen by the user)

		Deleted Files:
		Rplots.pdf
		depth.txt
		VCF_dict_tmp"""
echo $USAGE
}

pathToSampleDir=""
while getopts ":h:i:" opt; do 
	case $opt in
		h)
			usage
			exit 1
			;;
		i)
			pathToSampleDir=$OPTARG
			;;
		\?)
			echo "Invalid option: -$OPTARG" >&2
			exit 1
			;;
		:)
			echo "Option -$OPTARG requires an argument." >&2
			exit 1
			;;
	esac
done

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
myMtoolbox_out="${pathToSampleDir}/myMtoolbox_out"
mkdir $myMtoolbox_out
mv *.rg.* $myMtoolbox_out 
echo "/myMtoolbox_out holds the bam/sam outputs of the indel realignment and marking of duplicates"

rm Rplots.pdf
rm depth.txt
rm VCF_dict_tmp

echo ""
results="${pathToSampleDir}/results"
mkdir $results
mv *txt $results
mv *pdf $results
mv *vcf $results
echo "/myMtoolbox_out holds the bam/sam outputs of the indel realignment and marking of duplicates"


