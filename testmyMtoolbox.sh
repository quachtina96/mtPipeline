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
fasta_path=${mtoolbox_folder}data/ #might be something wrong here
mtdb_fasta=chr${ref}.fa 
hg19_fasta=/gpfs/group/stsi/genomes/GATK_bundle/hg19/ucsc.hg19.fasta 
samtoolsexe=/gpfs/group/stsi/methods/variant_calling/bwa_GATK/bin/samtools

cd /gpfs/home/quacht/test_myMtoolbox_ID18Father/tryMarkedBam

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

