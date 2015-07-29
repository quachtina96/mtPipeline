#!/bin/bash
set -e
set -o pipefail

cd /gpfs/home/quacht/toolbox2.0
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


for i in $(ls *rg.ra.marked.bam); do
echo "Converting ${i} to sam file..."
markedSam="$(echo ${i} | sed 's/.bam/.sam/')"
echo "this is the output: ${markedSam}"
java -Xmx4g \
-Djava.io.tmpdir=`pwd`/tmp \
-jar ${externaltoolsfolder}SamFormatConverter.jar \
INPUT="${i}" \
OUTPUT="$markedSam" \
TMP_DIR="`pwd`/tmp";
#commented out the following command (originally from MtoolBox, since it gave me errors during future bam analysis)
grep -v "^@" *marked.sam > "$(echo ${i} | sed 's/.sam/.norg.sam/')"
echo ${i} | sed 's/.sam/.norg.sam/'
done

#commented out the follwing section because I wanted to keep viewing the files while testing.
#mkdir MarkTmp
#mv *.ra.bam MarkTmp
#mv *.ra.marked.bam MarkTmp
#mv *marked.bam.marked.sam MarkTmp
#tar -czf MarkTmp.tar.gz MarkTmp
#rm -R MarkTmp

# ASSEMBLE CONTIGS, GET MT-TABLES...
echo ""
echo "##### ASSEMBLING MT GENOMES WITH ASSEMBLEMTGENOME..."
echo ""
echo "WARNING: values of tail < 5 are deprecated and will be replaced with 5"
echo ""	
#for each directory labeled as an output, 
for i in $(ls *rg.ra.marked.bam); do 
outhandle=$(echo ${i} | sed 's/.rg.ra.marked.sam//g')-mtDNAassembly; 
echo $outhandle
python /gpfs/home/quacht/toolbox/myAssembleMTgenome.py \
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
python VCFoutput.py -r ${ref}
