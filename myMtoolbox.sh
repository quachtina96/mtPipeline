#!/bin/bash
set -e
set -o pipefail

#need to test with the following line commented out
#source /gpfs/home/quacht/scripts/parameters.sh

usage()
{
	USAGE="""
	myMtoolbox.sh is a script that incorporates code from MToolBox's mitchondrial genome analysis pipeline for variant calling. 
	Written by Tina Quach, incorporating code from MToolBox (written by Domenico Simone, Claudia Calabrese and Maria Angela Diroma 2013-2014).
	
	Mandatory input options:

		-i	path to sample directory (should contain a mtExtract.csort.bam (a coordinate-sorted bam file of the extract mitchondrial region of an exome) 
			that results from running simplepipeline.py on the sample directory)

	Help options:

		-h	view usage

	"""
	echo "$USAGE"
}

sampleDir=""

while getopts ":h:i:" opt; do 
	case $opt in
		h)
			usage
			exit 1
			;;
		i)
			sampleDir=$OPTARG
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

cd $sampleDir
echo "Currently working in $sampleDir"

echo ""
echo "##### REALIGNING KNOWN INDELS WITH GATK INDELREALIGNER..."
echo ""

for i in *_mtExtractremap.csort.bam ; do 
echo "Currently working with ${i}..."

#First, add read groups to the file so IndelRealigner will be able to process the file
sampleName="${i//_exome_mtExtractremap.csort.bam/}"
echo "Adding read groups to the bam files";
java -Xmx2g \
-Djava.io.tmpdir=$(pwd)/tmp \
-jar ${picard}AddOrReplaceReadGroups.jar \
INPUT="${i}" \
OUTPUT="${sampleName}.rg.bam" \
SORT_ORDER=coordinate \
RGID="$sampleName" \
RGLB="$sampleName" \
RGPL="ILLUMINA" \
RGPU="L00" \
RGSM="${sampleName}_L00" ;

echo "..."
samtools index "${sampleName}.rg.bam"

echo "Realigning known indels for file ${i} using ${mtoolbox_folder}data/MITOMAP_HMTDB_known_indels.vcf as reference..."
java -Xmx4g \
-Djava.io.tmpdir=`pwd`/tmp \
-jar "$GATK" \
-T IndelRealigner \
-R ${mtoolbox_folder}/data/chr${ref}.fa \
-I "${sampleName}.rg.bam" \
-o "${sampleName}.rg.ra.bam" \
-targetIntervals ${mtoolbox_folder}data/intervals_file_RCRS.list \
-known ${mtoolbox_folder}data/MITOMAP_HMTDB_known_indels_RCRS.vcf \
-compress 0;

echo ""
echo "##### ELIMINATING PCR DUPLICATES WITH PICARDTOOLS MARKDUPLICATES..."
echo "" 
java -Xmx4g \
-Djava.io.tmpdir="$(pwd)/tmp" \
-jar "${picard}MarkDuplicates.jar" \
INPUT="${sampleName}.rg.ra.bam" \
OUTPUT="${sampleName}.rg.ra.marked.bam" \
METRICS_FILE="${sampleName}-metrics.txt" \
ASSUME_SORTED=true \
REMOVE_DUPLICATES=true \
TMP_DIR=$(pwd)/tmp; done >> ${pathToSampleDirs}log.txt

#Convert the marked.bam file to sam file for later processing (in myAssembleMTgenome.py)
for i in *marked.bam; do
samtools view -h $i > ${i}.sam
done

#remove read groups from the sam file (that has gone through adding of RG, indel realignment, 
#and marking of duplicates) for downstream analysis in myAssembleMTgenome.py
for i in *marked.bam.sam; do
grep -v "^@" $i > "${i//.bam.sam/.norg.sam}"
echo "${i//.bam.sam/.norg.sam}"
done

# ASSEMBLE CONTIGS, GET MT-TABLES...
echo ""
echo "##### ASSEMBLING MT GENOMES WITH MYASSEMBLEMTGENOME..."
echo ""
echo "WARNING: values of tail < 5 are deprecated and will be replaced with 5"
echo ""	
#for each directory labeled as an output, 
for i in *rg.ra.marked.bam; do 
outhandle=$(echo ${i} | sed 's/.rg.ra.marked.bam//g')-mtDNAassembly; 
echo $outhandle
python ${mtPipelineScripts}myAssembleMTgenome.py \
-i ${i} \
-o ${outhandle} \
-r ${fasta_path} \
-f ${mtdb_fasta} \
-a ${hg19_fasta} \
-s ${samtoolsexe} \
-FCP #${assembleMTgenome_OPTS}
done > logassemble.txt

echo ""
echo "##### GENERATING VCF OUTPUT #############"
python ${mtPipelineScripts}VCFoutput.py -s $sampleName >> ${pathToSampleDirs}log.txt

echo "VCF for $sampleName generated"



