#!/bin/bash

# This is the main file

set -e
set -o pipefail

./bashgatk.sh proc ./confs/kapa_hyperexome.sh ../../projects/#
echo "$(date): Processing completed" >> ./bashgatk.log
./bashgatk.sh somaSNP ./confs/kapa_hyperexome.sh ../../projects/#
echo "$(date): SomaticSNP ready" >> ./bashgatk.log
./bashgatk.sh proc ./confs/kapa_hyperexome.sh ../../projects/#
echo "$(date): Processing completed" >> ./bashgatk.log
./bashgatk.sh somaSNP ./confs/kapa_hyperexome.sh ../../projects/#
echo "$(date): SomaticSNP ready" >> ./bashgatk.log
#./bashgatk.sh anno $1 $2
#echo "$(date): Annotated" >> ./bashgatk.log
#./bashgatk.sh 2csv $1 $2
#echo "$(date): CSV created" >> ./bashgatk.log
