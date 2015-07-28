#!/bin/bash
#PBS -l walltime=2:00:00

#TODO: 
#define mtoolbox)folder (path to folder) e.g mtoolbox_folder=$(which $me | sed "s/$me//g")
#externaltoolsfolder

set -e
set -o pipefail

cd /gpfs/home/quacht/ID18exome/merge/chrRCRS_remap/
proj_path=`pwd`
#r -p $proj_path/results
#mkdir -p $proj_path/oldparts

module load samtools
module load python
module load bwa
module load bedtools
module load R

mainFileOut="$proj_path/pipelineOut_newRef.txt"
remapOut="$proj_path/remapOut_newRef.txt"

date
echo "Merging part.bam files into single exome.bam file for each sample..."
python /gpfs/home/quacht/ID18exome/merge/chrRCRS_remap/Pipe1-Merge.py > $mainFileOut
#mv *part.bam $proj_path/oldparts
echo "Calculating and analyzing the depth of coverage for the chrM region of the exome.bam for each sample"
echo "and extracting chrM region..."
python /gpfs/home/quacht/ID18exome/merge/chrRCRS_remap/Pipe2-getDepthAndExtract.py >> $mainFileOut
echo "Remapping the extracted chrM to rCRS..."
echo "See $remapOut for output"
python /gpfs/home/quacht/ID18exome/merge/chrRCRS_remap/Pipe3-RemapToRCRS.py > $remapOut
echo "Recalculating and analyze coverage for the sorted, remapped extracted chrM. "
python /gpfs/home/quacht/ID18exome/merge/chrRCRS_remap/Pipe4-getDepthsOfSortedRemapped.py >> $mainFileOut
echo "Check $mainFileOut for file ouput"
date

