#!/bin/bash

set -e
set -o pipefail

./bashgatk.sh procAmp $1 $2
echo "$(date): Processing completed" >> ./bashgatk.log
./bashgatk.sh deep $1 $2
echo "$(date): SomaticSNP ready" >> ./bashgatk.log
./bashgatk.sh anno $1 $2
echo "$(date): Annotated" >> ./bashgatk.log
./bashgatk.sh 2csv $1 $2
echo "$(date): CSV created" >> ./bashgatk.log
