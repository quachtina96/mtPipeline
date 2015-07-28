#!/bin/bash
#PBS -l walltime=2:00:00

#TODO: 
#define mtoolbox)folder (path to folder) e.g mtoolbox_folder=$(which $me | sed "s/$me//g")
#externaltoolsfolder

set -e
set -o pipefail

cd /gpfs/home/quacht/ID18exome/merge/
proj_path=`pwd`
mkdir -p $proj_path/results
mkdir -p $proj_path/oldparts

module load samtools
module load python
module load bwa
module load bedtools
module load R

fileoutput="$proj_path/pipelineOut.txt"

echo "Merging part.bam files into single exome.bam file for each sample..."
python /gpfs/home/quacht/ID18exome/merge/Pipe1-Merge.py > $fileoutput
#mv *part.bam $proj_path/oldparts
echo "Calculating and analyzing the depth of coverage for the chrM region of the exome.bam for each sample"
echo "and extracting chrM region..."
python /gpfs/home/quacht/ID18exome/merge/Pipe2-getDepthAndExtract.py >> $fileoutput
echo "Remapping the extracted chrM to rCRS..."
echo "See remapOut.txt for output"
python /gpfs/home/quacht/ID18exome/merge/Pipe3-Remap.py > remapOut.txt
echo "Recalculating and analyze coverage for the sorted, remapped extracted chrM. "
python /gpfs/home/quacht/ID18exome/merge/Pipe4-getDepthsOfSortedRemapped.py >> $fileoutput


