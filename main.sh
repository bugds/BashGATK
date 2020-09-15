#!/bin/bash

set -e
set -o pipefail

./bashgatk.sh proc $1
echo "$(date): Processing completed" >> ./bashgatk.log
./bashgatk.sh somaSNP $1
echo "$(date): SomaticSNP ready" >> ./bashgatk.log
./bashgatk.sh anno $1
echo "$(date): Annotated" >> ./bashgatk.log
./bashgatk.sh 2csv $1
echo "$(date): CSV created" >> ./bashgatk.log
