#!/bin/bash
set -e
set -o pipefail


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
fasta_path=/gpfs/home/quacht/data/ 
mtdb_fasta=chrRCRS.fa 
hg19_fasta=/gpfs/group/stsi/genomes/GATK_bundle/hg19/ucsc.hg19.fasta 
samtoolsexe=/gpfs/group/stsi/methods/variant_calling/bwa_GATK/bin/samtools

cd /gpfs/home/quacht/toolbox1.0
currentDir=$(pwd)
echo $currentDir

echo ""
echo "##### GENERATING VCF OUTPUT..."
# ... AND VCF OUTPUT
python /gpfs/home/quacht/scripts/VCFoutput.py -r ${ref}

