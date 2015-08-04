#!/bin/bash
set -e
set -o pipefail


cd /gpfs/home/quacht/ID18
pwd

for sampleDir in *; do
if [ "$sampleDir" != "ID18_Father" ] &&  [ "$sampleDir" != "VCF" ] && [ "$sampleDir" != "log.txt" ]; then
echo $sampleDir
fi
done