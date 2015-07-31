#!/bin/bash
set -e
set -o pipefail

#NOTE: THIS HAS BEEN WRITTEN AS IF I WILL HAVE EVERYTHING RELEVANT PACKAGED INTO A MTPIPELINE FOLDER LIKE MTOOLBOX
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
 
		-i	path to input folder.
		-o	path to output folder.
		-m  path to mtPipeline folder.
		-a	options for assembleMTgenome script [see assembleMTgenome.py -h for details]
		-c	options for mt-classifier script [see mt-classifier.py -h for details]

	Help options:

		-h	show this help message

	"""
	echo "$USAGE"
}

source /gpfs/home/quacht/scripts/parameters.sh
#statements

export mtPipelineScripts=/gpfs/home/quacht/scripts/
export mtoolbox_folder=/gpfs/home/quacht/MToolBox/
export externaltoolsfolder=/gpfs/home/quacht/MToolBox/ext_tools/
export ref="RCRS"
export fasta_path=${mtoolbox_folder}data/ #might be something wrong here
export mtdb_fasta=chr${ref}.fa 
export hg19_fasta=/gpfs/group/stsi/genomes/GATK_bundle/hg19/ucsc.hg19.fasta 
export samtoolsexe=/gpfs/group/stsi/methods/variant_calling/bwa_GATK/bin/samtools


while getopts ":h:a:c:i:o:m:" opt; do 
	case $opt in
		h)
			usage
			exit 1
			;;
		a)
			assembleMTgenome_OPTS=$OPTARG
			;;
		c)
			mt_classifier_OPTS=$OPTARG
			;;
		i)
			input_path=$OPTARG
			;;
		o)
			output_name=$OPTARG
			;;	
		m)
			mtPipeFolder=$OPTARG	
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

cd "${input_path}"
#python simplepipe.py -m ${mtPipeFolder} -i ${input_path} > log-simplepipe.txt
python ${mtPipelineScripts}simplepipe.py -i ${input_path} > log-simplepipe.txt
bash ${mtPipelineScripts}myMtoolbox.sh -i  ${input_path} > log-myMtoolbox.txt

