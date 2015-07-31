#!/bin/bash
set -e
set -o pipefail

source parameters.sh

export mtPipeFolder=/gpfs/home/quacht/mtPipeline
export mtoolbox_folder=/gpfs/home/quacht/MToolBox/
export externaltoolsfolder=/gpfs/home/quacht/MToolBox/ext_tools/
export ref="RCRS"
export fasta_path=${mtoolbox_folder}data/ #might be something wrong here
export mtdb_fasta=chr${ref}.fa 
export hg19_fasta=/gpfs/group/stsi/genomes/GATK_bundle/hg19/ucsc.hg19.fasta 
export samtoolsexe=/gpfs/group/stsi/methods/variant_calling/bwa_GATK/bin/samtools

usage()
{
	USAGE="""
	myMtoolbox.sh is a script that incorporates code from MToolBox's mitchondrial genome analysis pipeline for variant calling. 
	Written by Tina Quach, incorporating code from MToolBox (written by Domenico Simone, Claudia Calabrese and Maria Angela Diroma 2013-2014).
	
	Mandatory input options:

		-i	path to input folder.
		-m  path to mtPipeline folder.

	Help options:

		-h	view usage

	"""
	echo "$USAGE"
}

sampleDir=""

while getopts ":h:i:m" opt; do 
	case $opt in
		h)
			usage
			exit 1
			;;
		i)
			sampleDir=$OPTARG
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

cd $sampleDir
echo "Currently working in $sampleDir"

echo ""
echo "##### REALIGNING KNOWN INDELS WITH GATK INDELREALIGNER..."
echo ""
for i in $(ls *_mtExtractremap.csort.bam); do 
echo "Currently working with ${i}..."
sampleName="$(echo ${i} | sed 's/_exome_mtExtractremap.csort.bam//')"
echo "Adding read groups to the bam files";
java -Xmx2g \
-Djava.io.tmpdir=`pwd`/tmp \
-jar /opt/applications/picard/current/AddOrReplaceReadGroups.jar \
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

#java -jar /opt/applications/picard/current/CreateSequenceDictionary.jar R=/gpfs/home/quacht/ID18exome/rCRS.fasta O=/gpfs/home/quacht/ID18exome/rCRS.dict

echo "Realigning known indels for file ${i} using ${mtoolbox_folder} data/MITOMAP_HMTDB_known_indels.vcf as reference..."
java -Xmx4g \
-Djava.io.tmpdir=`pwd`/tmp \
-jar /opt/applications/gatk/3.3-0/GenomeAnalysisTK.jar \
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
-jar /opt/applications/picard/current/MarkDuplicates.jar \
INPUT="${sampleName}.rg.ra.bam" \
OUTPUT="${sampleName}.rg.ra.marked.bam" \
METRICS_FILE="${sampleName}-metrics.txt" \
ASSUME_SORTED=true \
REMOVE_DUPLICATES=true \
TMP_DIR=`pwd`/tmp; done

for i in $(ls *marked.bam); do
samtools view -h $i > ${i}.sam
done

for i in $(ls *marked.bam.sam); do

grep -v "^@" $i > "$(echo ${i} | sed 's/.bam.sam/.norg.sam/')"
echo ${i} | sed 's/.bam.sam/.norg.sam/'
done

# ASSEMBLE CONTIGS, GET MT-TABLES...
echo ""
echo "##### ASSEMBLING MT GENOMES WITH ASSEMBLEMTGENOME..."
echo ""
echo "WARNING: values of tail < 5 are deprecated and will be replaced with 5"
echo ""	
#for each directory labeled as an output, 
for i in $(ls *rg.ra.marked.bam); do 
outhandle=$(echo ${i} | sed 's/.rg.ra.marked.bam//g')-mtDNAassembly; 
echo $outhandle
python /gpfs/home/quacht/scripts/myAssembleMTgenome.py \
-i ${i} \
-o ${outhandle} \
-r ${fasta_path} \
-f ${mtdb_fasta} \
-a ${hg19_fasta} \
-s ${samtoolsexe} \
-FCP #${assembleMTgenome_OPTS}
done > logassemble.txt

echo ""
echo "##### GENERATING VCF OUTPUT..."
# ... AND VCF OUTPUT
python /gpfs/home/quacht/scripts/VCFoutput.py -r ${ref}


