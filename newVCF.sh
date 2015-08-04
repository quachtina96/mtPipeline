#!/bin/bash
set -e
set -o pipefail

#!/bin/bash
set -e
set -o pipefail

source /gpfs/home/quacht/scripts/parameters.sh

cd /gpfs/home/quacht/testVCFmerge

for folder in *; do
bash ${mtPipelineScripts}myMtoolbox.sh -i $folder > log.txt
done 