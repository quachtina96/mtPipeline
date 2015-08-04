#!/bin/bash

# executables, bwa5=bwa-0.5.8c, bwa6=bwa-0.6.2, bwa7=bwa-0.7.4
BIN=/gpfs/group/stsi/methods/variant_calling/bwa_GATK/bin
BWA5=$BIN/bwa5
BWA6=$BIN/bwa6
BWA7=$BIN/bwa7
module load bwa
samtoolsexe=$BIN/samtools
module load samtools
#GATK 3.3.0
GATK=/opt/applications/gatk/3.3-0/GenomeAnalysisTK.jar
JAVA=/usr/bin/java
module load java/1.7.0_21
module load R
module load bedtools
module load python
picard=/opt/applications/picard/current/

#MToolBox bundle
mtoolbox_folder=/gpfs/home/quacht/MToolBox/
externaltoolsfolder=/gpfs/home/quacht/MToolBox/ext_tools/
fasta_path=${mtoolbox_folder}data/
ref="RCRS"
mtdb_fasta=chr${ref}.fa 
hg19_fasta=/gpfs/group/stsi/genomes/GATK_bundle/hg19/ucsc.hg19.fasta 

#mtPipeline
mtPipelineFolder=/gpfs/home/quacht/
mtPipelineScripts=${mtPipelineFolder}scripts/


