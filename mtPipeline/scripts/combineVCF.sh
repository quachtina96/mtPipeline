#!/bin/bash
set -e
set -o pipefail

usage()
{
	USAGE="""
	combineVCF.sh is a script that will combine VCF files. It is designed to combine 3 samples.

		Mandatory input options:

		-i	path to VCFDir.
		-p  path to parameters.sh
	Help options:

		-h	view usage

	"""
	echo "$USAGE"
}

while getopts ":h:i:p:" opt; do 
	case $opt in
		h)
			usage
			exit 1
			;;
		i)
			VCFDir=$OPTARG
			;;
		p) 
			parameters=$OPTARG
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

source "$parameters"
echo $VCFDir
cd $VCFDir
echo "Currently working in $VCFDir"

echo ""
echo "##### COMBINING VCF FILES..."
echo ""

VCFarray=(*VCF.vcf)
echo "${VCFarray[0]}"
echo "${VCFarray[1]}"
echo "${VCFarray[2]}"

java -jar /opt/applications/gatk/3.3-0/GenomeAnalysisTK.jar \
   -T CombineVariants \
   -R ${fasta_path}chrRCRS.fa \
   --variant:father "${VCFDir}/${VCFarray[0]}" \
   --variant:mother "${VCFDir}/${VCFarray[1]}" \
   --variant:proband "${VCFDir}/${VCFarray[2]}" \
	  -o combined.vcf \
