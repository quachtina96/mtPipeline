#!/bin/bash
set -e
set -o pipefail

#do i need to source parameter.sh here?
source /gpfs/home/quacht/scripts/parameters.sh

usage()
{
	USAGE="""
	myMtoolbox.sh is a script that incorporates code from MToolBox's mitchondrial genome analysis pipeline for variant calling. 
	Written by Tina Quach, incorporating code from MToolBox (written by Domenico Simone, Claudia Calabrese and Maria Angela Diroma 2013-2014).
	
	Mandatory input options:

		-i	path to VCFDir.

	Help options:

		-h	view usage

	"""
	echo "$USAGE"
}

while getopts ":h:i:m" opt; do 
	case $opt in
		h)
			usage
			exit 1
			;;
		i)
			VCFDir=$OPTARG
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

cd $VCFDir
echo "Currently working in $VCFDir"

echo ""
echo "##### COMBINING VCF FILES..."
echo ""

VCFarray=(*vcf)
echo "${VCFarray[0]}"
echo "${VCFarray[1]}"
echo "${VCFarray[2]}"

java -jar GenomeAnalysisTK.jar \
   -T CombineVariants \
   -R chrRCRS.fa \
   --variant "${VCFarray[0]}" \
   --variant "${VCFarray[1]}" \
   --variant "${VCFarray[2]}" \
      -o combined.vcf \
   -assumeIdentiticalSamples