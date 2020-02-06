#!/bin/bash

set -e
set -o pipefail

export inputFolder="$(realpath $1)/"
export outputFolder="$(realpath $2)/"
export gatk=/opt/gatk-4.1.4.1/gatk

export regions=/home/bioinfuser/NGS/Reference/intervals/2020_02_02/regions.bed
export refFasta=/home/bioinfuser/NGS/Reference/hg38/hg38.fasta
export refDict=/home/bioinfuser/NGS/Reference/hg38/hg38.dict
export gnomad=/home/bioinfuser/NGS/Reference/hg38/Mutect2/gnomad.exomes.r2.1.1.sites.liftover_grch38.vcf.bgz
export variantsForContamination=/home/bioinfuser/NGS/Reference/hg38/Mutect2/variants_for_contamination.vcf

export javaOpt="-Xms3000m"

function makeDirectory {
    local newDirectory=$1

    if !( test -d "${outputFolder}${newDirectory}/" ); then
        mkdir "${outputFolder}${newDirectory}/"
    fi
}

function mutect2 {
    local files="${outputFolder}recalibrated/*.bam"
    makeDirectory mutect2

    for bam in $files; do
        local bamName = $(basename -- ${bam})
        makeDirectory mutect2/${bamName}

        touch ${outputFolder}mutect2/${bamName}/bamout.bam
        touch ${outputFolder}mutect2/${bamName}/f1r2.tar.gz

        $gatk --java-options "${javaOpt}" \
          GetSampleName \
            -R refFasta \
            -I $bam \
            -O ${outputFolder}mutect2/${bamName}/${bamName}.txt \
            -encode

        $gatk --java-options "${javaOpt}" \
          Mutect2 \
            -R refFasta \
            -I $bam \
            -tumor 'cat ${outputFolder}mutect2/${bamName}/${bamName}.txt' \
            --germline-resource $gnomad \
            -L $regions \
            -O ${outputFolder}mutect2/${bamName}/${bamName}.unfiltered.vcf \
            --bam-output bamout.bam \
            --f1r2-tar-gz f1r2.tar.gz
        
        $gatk --java-options "${javaOpt}" \
          GetPileupSummaries \
            -R refFasta \
            -I $bam \
            --interval-set-rule INTERSECTION -L $regions \
            -V $variantsForContamination \
            -L $variantsForContamination \
            -O ${outputFolder}mutect2/${bamName}/${bamName}.pileups.table
    done
}

# Merging VCFs and .bams is for WGS, when it's beneficial
# to use intervals, therefore paralleling the process

function mergeVcfs {
    local files=${outputFolder}mutect2/*/*.unfiltered.vcf
    local inputFiles=$(printf -- "--INPUT %s " $files)

    $gatk --java-options "${javaOpt}" \
      MergeVcfs \
        $inputFiles \
        -O ${outputFolder}mutect2/unfiltered.vcf
}

function mergeBamOuts {
    local files=${outputFolder}mutect2/*/bamout.bam
    local inputFiles=$(printf -- "--INPUT %s " $files)

    $gatk --java-options "${javaOpt}" \
      GatherBamFiles \
        $inputFiles \
        -O ${outputFolder}mutect2/unsorted.out.bam
}

function mergePileupSummaries {
    local files=${outputFolder}mutect2/*/*.pileups.table
    local inputFiles=$(printf -- "--INPUT %s " $files)

    $gatk --java-options "${javaOpt}" \
      GatherPileupSummaries \
        --sequence-dictionary $refDict
        $inputFiles
        -O ${outputFolder}mutect2/gathered.pileups.table

}

function learnReadOrientationModel {
    local files=${outputFolder}mutect2/*/f1r2.tar.gz
    local inputFiles=$(printf -- "--INPUT %s " $files)

    $gatk --java-options "${javaOpt}" \
      LearnReadOrientationModel \
        $inputFiles
        -O artifact-priors.tar.gz
}

function calculateContamination {
    $gatk --java-options "${javaOpt}" \
      CalculateContamination 
        -I ${outputFolder}mutect2/gathered.pileups.table \
        -O contamination.table \
        --tumor-segmentation ${outputFolder}mutect2/segments.table
}

mutect2

# There is no need to run any of these functions for all files
# in the run:
# https://github.com/broadinstitute/gatk/tree/master/scripts/mutect2_wdl
# The repository above indicates that there are no difference
# in initiation of Mutect2 pipeline functions