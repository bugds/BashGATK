#!/bin/bash

set -e
set -o pipefail

export command=$1
export config=$2
export outputFolder="$(realpath $3)/"

source $config

if [ $command == 'proc' ]; then
    bash ${scriptsDirectory}/parallelProcessing.sh
elif [ $command == 'procAmp' ]; then
    bash ${scriptsDirectory}/parallelProcessingAmpliconBased.sh
elif [ $command == 'procWGS' ]; then
    bash ${scriptsDirectory}/processingWGS.sh
elif [ $command == 'somaSNP' ]; then
    bash ${scriptsDirectory}/parallelSomaticSNP.sh
elif [ $command == 'germSNP' ]; then
    bash ${scriptsDirectory}/parallelGermlineSNP.sh
elif [ $command == 'anno' ]; then
    bash ${scriptsDirectory}/annotation.sh
elif [ $command == 'qgen' ]; then
    bash ${scriptsDirectory}/annotation_qiagen.sh
elif [ $command == 'btil' ]; then
    bash ${scriptsDirectory}/bed_to_interval_list.sh
elif [ $command == 'cvfc' ]; then
    bash ${scriptsDirectory}/create_variants_for_contamination.sh
elif [ $command == '2csv' ]; then
    python3 ${scriptsDirectory}/goCsv.py $outputFolder
elif [ $command == '2csvq' ]; then
    python3 ${scriptsDirectory}/goCsvQiagen.py $outputFolder
elif [ $command == 'cnvk' ]; then
    bash ${scriptsDirectory}/run_cnvkit.sh
elif [ $command == 'deep' ]; then
    bash ${scriptsDirectory}/deepvariant.sh
elif [ $command == 'deep19' ]; then
    bash ${scriptsDirectory}/deepvariant.sh
elif [ $command == 'cnvk19' ]; then
    bash ${scriptsDirectory}/run_cnvkit.sh
elif [ $command == '2rus' ]; then
    python3 ${scriptsDirectory}/2rus.py $outputFolder
elif [ $command == 'cint' ]; then
    bash ${scriptsDirectory}/createIntervals.sh $outputFolder
elif [ $command == 'aved' ]; then
    bash ${scriptsDirectory}/averageDepth.sh $outputFolder
fi
