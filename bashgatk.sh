#!/bin/bash

set -e
set -o pipefail

export command=$1
export outputFolder="$(realpath $2)/"

source ./conf.sh

echo ${scriptsDirectory}

if [ $command == 'proc' ]; then
    bash ${scriptsDirectory}/parallelProcessing.sh
elif [ $command == 'procWGS' ]; then
    bash ${scriptsDirectory}/processingWGS.sh
elif [ $command == 'somaSNP' ]; then
    bash ${scriptsDirectory}/parallelSomaticSNP.sh
elif [ $command == 'germSNP' ]; then
    bash ${scriptsDirectory}/parallelGermlineSNP.sh
elif [ $command == 'anno' ]; then
    bash ${scriptsDirectory}/annotation.sh
elif [ $command == 'btil' ]; then
    bash ${scriptsDirectory}/bed_to_interval_list.sh
elif [ $command == 'cvfc' ]; then
    bash ${scriptsDirectory}/create_variants_for_contamination.sh
elif [ $command == '2csv' ]; then
    python3 ${scriptsDirectory}/goCsv.py $outputFolder
elif [ $command == 'cnvk' ]; then
    bash ${scriptsDirectory}/run_cnvkit.sh
elif [ $command == 'deep' ]; then
    bash ${scriptsDirectory}/deepvariant.sh
fi
