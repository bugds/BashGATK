#!/bin/bash

# This is the main file

set -e
set -o pipefail

conf='./confs/example.sh'

./bashgatk.sh kumi $conf
echo "$(date): Processing completed" >> ./bashgatk.log
./bashgatk.sh deep $conf
echo "$(date): DeepVariant ready" >> ./bashgatk.log
./bashgatk.sh somaSNP $conf
echo "$(date): SomaticSNP ready" >> ./bashgatk.log
./bashgatk.sh anno $conf
echo "$(date): Annotation completed" >> ./bashgatk.log
./bashgatk.sh 2csv $conf
echo "$(date): CSV created" >> ./bashgatk.log