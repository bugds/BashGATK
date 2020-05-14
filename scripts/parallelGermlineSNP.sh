#!/bin/bash 

set -e
set -o pipefail

export regions=/home/bioinfuser/NGS/Reference/intervals/2020_02_02/capture_targets.interval_list
export parallelJobs=5

export javaOpt="-Xms3000m"

function makeDirectory {
    local newDirectory=$1

    if !( test -d "${outputFolder}${newDirectory}/" ); then
        mkdir "${outputFolder}${newDirectory}/"
    fi
}

function parallelRun {  
    local files=$2

    export -f $1
    parallel -j $parallelJobs $1 ::: $files
}

function haplotypeCaller {
    local bamName=$(basename -- ${1} | cut -f 1 -d '.')

    $gatk --java-options "${javaOpt}" \
      HaplotypeCaller \
        -R $refFasta \
        -I $1 \
        -L $regions \
        -O ${outputFolder}haplotype_caller/${bamName}.vcf \
        -ERC GVCF \
        -contamination 0 \
        -G StandardAnnotation \
        -G AS_StandardAnnotation \
        -G StandardHCAnnotation
}

makeDirectory haplotype_caller
parallelRun haplotypeCaller "${outputFolder}recalibrated/*.bam"
sleep 1

