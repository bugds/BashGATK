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

export numCpu=16
export regions=/home/bioinfuser/NGS/Reference/intervals/2020_02_02/capture_targets.bed
export binVersion='0.10.0'
export refFolder=/home/bioinfuser/NGS/Reference/
export refFastaPathPart=hg38/hg38.fasta
export regionsPathPart=intervals/2020_02_02/cnvkit_targets.bed

function makeDirectory {
    local newDirectory=$1

    if !( test -d "${outputFolder}${newDirectory}/" ); then
        mkdir "${outputFolder}${newDirectory}/"
    fi
}

function runDeepvariant {
    makeDirectory "deepvariant"
    
    files=${outputFolder}recalibrated/*.bam

    for bam in $files; do
        local bamName=$(basename -- $bam | cut -f1 -d '.')
        local baseName=$(basename -- $bam)
        
        sudo docker run \
          -v ${outputFolder}recalibrated:"/input" \
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
