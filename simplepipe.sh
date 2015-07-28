#!/bin/bash
#PBS -l walltime=2:00:00

#TODO: 
#define mtoolbox)folder (path to folder) e.g mtoolbox_folder=$(which $me | sed "s/$me//g")
#externaltoolsfolder

set -e
set -o pipefail

cd /gpfs/home/quacht/scripts/
proj_path=`pwd`

module load samtools
module load python
module load bwa
module load bedtools
module load R

fileoutput="$proj_path/test_simplepipe_ID18Mother.txt"
python /gpfs/home/quacht/scripts/test_simplepipe_ID18Mother.py > $fileoutput 



