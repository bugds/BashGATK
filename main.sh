#!/bin/bash

set -e
set -o pipefail

./bashgatk.sh proc ./confs/conf_bykova_more.sh ../../projects/2022_03_Bykova
echo "$(date): Processing completed" >> ./bashgatk.log
./bashgatk.sh somaSNP ./confs/conf_bykova_more.sh ../../projects/2022_03_Bykova
echo "$(date): SomaticSNP ready" >> ./bashgatk.log
./bashgatk.sh proc ./confs/conf_bykova_less.sh ../../projects/2022_03_Bykova_less
echo "$(date): Processing completed" >> ./bashgatk.log
./bashgatk.sh somaSNP ./confs/conf_bykova_less.sh ../../projects/2022_03_Bykova_less
echo "$(date): SomaticSNP ready" >> ./bashgatk.log
#./bashgatk.sh anno $1 $2
#echo "$(date): Annotated" >> ./bashgatk.log
#./bashgatk.sh 2csv $1 $2
#echo "$(date): CSV created" >> ./bashgatk.log
