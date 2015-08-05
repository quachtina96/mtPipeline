#!/bin/bash
set -e
set -o pipefail

usage()
{
	USAGE="""
	mtPipeline: a pipeline for heteroplasmy annotation and accurate functional analysis of mitochondrial variants from high throughput sequencing data.
	Written by Tina Quach, incorporating code from MToolBox (written by Domenico Simone, Claudia Calabrese and Maria Angela Diroma 2013-2014).
	
	You must run the mtPipeline command on only on bam file types.
	The input folder should contain a folder for each sample. 
	Each sample folder is expected to contain one or more bam files for the particular sample. Those bam files will be merged to form a single bam file 
	for the sample.

	uses only rCRS

	Input & workflow execution options (must include -i and -m):
 
		-i	path to input folder (containing sample directories). INCLUDE the last backslash.
		-p  path to parameters.sh
		-a	options for assembleMTgenome script [see assembleMTgenome.py -h for details]

	Help options:

		-h	show this help message

	"""
	echo "$USAGE"
}


while getopts ":h:a:c:i:o:p:" opt; do 
	case $opt in
		h)
			usage
			exit 1
			;;
		a)
			assembleMTgenome_OPTS=$OPTARG
			;;
		i)
			pathToSampleDirs=$OPTARG
			;;
		p)
			pathToParameters=$OPTARG
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


source $pathToParameters

echo ""
echo "Check python version... (2.7 required)"
min=$(python -c "import sys; print (sys.version_info[:])[1]")
maj=$(python -c "import sys; print (sys.version_info[:])[0]")
if [[ $maj != 2 ]] || [[ $min != 7 ]]
then
echo "You need Python2.7 in order to run MToolBox. Abort."
exit 1
else
echo "OK."
echo ""
fi

cd "$pathToSampleDirs"
pwd

for sampleDir in *; do
if  [ "$sampleDir" != "VCF" ] && [ "$sampleDir" != "log.txt" ]; then
#run the analysis on single sample
echo "Working with $sampleDir"
pathToSampleDir="${pathToSampleDirs}${sampleDir}"

#run simplepipe.py in order to merge the part.bam files, analyze coverage, extract chrM, remap the chrM to rCRS, and recalculate coverage 
python ${mtPipelineScripts}simplepipe.py -i ${pathToSampleDir} -m ${mtPipelineFolder} >> ${pathToSampleDirs}log.txt
#run myMtoolbox.sh to further process tge bam file resulting from above (Add RG, indel realign, mark duplicates, assemble MTgenome, variant call)
bash ${mtPipelineScripts}myMtoolbox.sh -i  ${pathToSampleDir} -p ${pathToParameters}>> ${pathToSampleDirs}log.txt



#clean up the sample directory
echo "############ ORGANIZING OUTPUTS ##############"
bash "${mtPipelineScripts}cleanUp.sh" -i "${pathToSampleDir}" >> ${pathToSampleDirs}log.txt
pwd

#copy the vcf for that sample to a VCF folder for the sample set for future VCF analysis
echo "############### COPY VCF ###############"
cd $pathToSampleDirs
pwd
if [ ! -d "${pathToSampleDirs}VCF" ] ; then
	 mkdir "${pathToSampleDirs}VCF"
fi

cd ${pathToSampleDir}/results/
cp *.vcf "${pathToSampleDirs}VCF"

fi

done

#combine the VCFs for the cohort
cd $pathToSampleDirs
bash ${mtPipelineScripts}combineVCF.sh -i "${pathToSampleDirs}VCF" -p "$pathToParameters" >> log.txt