#!/bin/bash

set -e
set -o pipefail

# ./bashgatk.sh proc $1
# ./bashgatk.sh somaSNP $1
./bashgatk.sh anno $1
# ./bashgatk.sh 2csv $1
