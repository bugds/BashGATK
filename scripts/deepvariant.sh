#!/bin/bash

set -e
set -o pipefail

# in this script:
# 1) open docker container in current wd which has standard directory structure
# 2) make an automatic script inside this container which runs by itself,
#    and processes all bam files:
#    https://stackoverflow.com/questions/47671589/how-can-i-run-script-automatically-after-docker-container-startup
# 3) after container is running, probably no commands could be driven from the
#    initial bash terminal

# Maybe get pisces, and get results of these 3 running through SURVIVOR of some kind

function makeDirectory {
    local newDirectory=$1

    if !( test -d "${outputFolder}${newDirectory}/" ); then
        mkdir "${outputFolder}${newDirectory}/"
    fi
}

function runDeepvariant {
    makeDirectory "deepvariant"
    
    files=${outputFolder}sorted/*.bam

    for bam in $files; do
        local bamName=$(basename -- $bam | cut -f1 -d '.')
        local baseName=$(basename -- $bam)
        
        docker run \
          -v ${outputFolder}sorted:"/input" \
          -v ${outputFolder}deepvariant:"/output" \
          -v ${refFolder}:"/reference" \
          google/deepvariant:${binVersion} \
          /opt/deepvariant/bin/run_deepvariant \
          --model_type=WES \
          --ref=/reference/${refFastaPathPart} \
          --reads=/input/${baseName} \
          --regions /reference/${regionsPathPart} \
          --output_vcf=/output/${bamName}.vcf.gz \
          --output_gvcf=/output/${bamName}.g.vcf.gz \
          --num_shards=${numCpu}
    done
}

runDeepvariant
