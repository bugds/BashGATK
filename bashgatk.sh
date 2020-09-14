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
<<<<<<< HEAD
elif [ $command == 'germSNP' ]; then
    bash ${scriptsDirectory}/parallelGermlineSNP.sh
elif [ $command == 'anno' ]; then
    bash ${scriptsDirectory}/annotation.sh
elif [ $command == 'btil' ]; then
=======
elif [ $cmd == 'germSNP' ]; then
    echo 'germSNP not implemented!'
elif [ $cmd == 'deep' ]; then
    bash ${scriptsDirectory}/deepvariant.sh
elif [ $cmd == 'anno' ]; then
    bash ${scriptsDirectory}/annotation.sh
elif [ $cmd == '2csv' ]; then
    python3 ${scriptsDirectory}/goCsv.py $outputFolder
elif [ $cmd == 'btil' ]; then
>>>>>>> 9e996baaa2d57bfaad090b553cb8fe9a6ae84bea
    bash ${scriptsDirectory}/bed_to_interval_list.sh
elif [ $command == 'cvfc' ]; then
    bash ${scriptsDirectory}/create_variants_for_contamination.sh
<<<<<<< HEAD
elif [ $command == '2csv' ]; then
    python3 ${scriptsDirectory}/goCsv.py $outputFolder
elif [ $command == 'cnvk' ]; then
    bash ${scriptsDirectory}/run_cnvkit.sh
elif [ $command == 'deep' ]; then
    bash ${scriptsDirectory}/deepvariant.sh
=======
elif [ $cmd == 'cnvk' ]; then
    bash ${scriptsDirectory}/run_cnvkit.sh
>>>>>>> 9e996baaa2d57bfaad090b553cb8fe9a6ae84bea
fi
