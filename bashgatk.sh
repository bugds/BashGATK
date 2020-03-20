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
    echo '2'
elif [ $cmd == 'somaSNP' ]; then
    bash ${scriptsDirectory}/parallelSomaticSNP.sh
elif [ $cmd == 'anno' ]; then
    bash ${scriptsDirectory}/annotation.sh
elif [ $cmd == 'btil' ]; then
    echo '5'
elif [ $cmd == 'cvfc' ]; then
    echo '6'
elif [ $cmd == '2csv' ]; then
    python3 ${scriptsDirectory}/goCsv.py $outputFolder
fi
