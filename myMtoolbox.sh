#!/bin/bash
set -e
set -o pipefail

cd /gpfs/home/quacht/checkMarkedDup
module load samtools
module load python
module load bwa
module load bedtools
module load R
module load java/1.7.0_21  
	#statements
mtoolbox_folder=/gpfs/home/quacht/MToolBox/
externaltoolsfolder=/gpfs/home/quacht/MToolBox/ext_tools/
ref="RCRS"
fasta_path=${mtoolbox_folder}data/ #might be something wrong here
mtdb_fasta=chr${ref}.fa 
hg19_fasta=/gpfs/group/stsi/genomes/GATK_bundle/hg19/ucsc.hg19.fasta 
samtoolsexe=/gpfs/group/stsi/methods/variant_calling/bwa_GATK/bin/samtools

#pwd and cd to correct folder
echo ""
echo "##### REALIGNING KNOWN INDELS WITH GATK INDELREALIGNER..."
echo ""
for i in $(ls *ID18_Father_exome_mtExtractremap.csort.bam); do 
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
-compress 0;done

for i in $(ls *bam);
do samtools view -h $i > ${i}.sam
done
