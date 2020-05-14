#!/bin/bash

set -e
set -o pipefail

cmd=$1
export outputFolder="$(realpath $2)/"
export scriptsDirectory=/home/bioinfuser/NGS/Pipelines/scripts

export gatk=/opt/gatk-4.1.5.0/gatk
export samtools=/opt/gatk4-data-processing/samtools-1.3.1/samtools
export picard=/opt/gatk4-data-processing/picard-2.16.0/picard.jar
export bwa=/opt/gatk4-data-processing/bwa-0.7.15/bwa

export refFasta=/home/bioinfuser/NGS/Reference/hg38/hg38.fasta
export refDict=/home/bioinfuser/NGS/Reference/hg38/hg38.dict

if [ $cmd == 'proc' ]; then
    bash ${scriptsDirectory}/parallelProcessing.sh
elif [ $cmd == 'procWGS' ]; then
    echo 'procWGS not implemented!'
elif [ $cmd == 'somaSNP' ]; then
    bash ${scriptsDirectory}/parallelSomaticSNP.sh
elif [ $cmd == 'germSNP' ]; then
    bash ${scriptsDirectory}/parallelGermlineSNP.sh
elif [ $cmd == 'anno' ]; then
    bash ${scriptsDirectory}/annotation.sh
elif [ $cmd == 'btil' ]; then
    bash ${scriptsDirectory}/bed_to_interval_list.sh
elif [ $cmd == 'cvfc' ]; then
    bash ${scriptsDirectory}/create_variants_for_contamination.sh
elif [ $cmd == '2csv' ]; then
    python3 ${scriptsDirectory}/goCsv.py $outputFolder
elif [ $cmd == 'cnvk' ]; then
    bash ${scriptsDirectory}/run_cnvkit.sh
elif [ $cmd == 'deep' ]; then
    bash ${scriptsDirectory}/deepvariant.sh
fi
