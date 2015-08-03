#!/bin/bash

set -e
set -o pipefail


cd /gpfs/home/quacht/testVCFmerge
pwd

VCFarray=(*vcf)
echo $VCFarray
echo "${VCFarray[0]}"
echo "${VCFarray[1]}"
echo "${VCFarray[2]}"