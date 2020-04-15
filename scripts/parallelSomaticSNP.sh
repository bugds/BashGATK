#!/bin/bash

# Multisample Mutect2 usage described here: 
# https://github.com/broadinstitute/gatk/blob/master/scripts/mutect2_wdl/mutect2_multi_sample.wdl

set -e
set -o pipefail

export regions=/home/bioinfuser/NGS/Reference/intervals/2020_02_02/capture_targets.interval_list
export refImg=/home/bioinfuser/NGS/Reference/hg38/hg38.fasta.img
export gnomad=/home/bioinfuser/NGS/Reference/intervals/2020_02_02/additional/AFonly.vcf
export variantsForContamination=/home/bioinfuser/NGS/Reference/intervals/2020_02_02/additional/variants_for_contamination.vcf
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

function mutect2 {
    local bamName=$(basename -- ${1} | cut -f 1 -d '.')
    makeDirectory mutect2/${bamName}

    touch ${outputFolder}mutect2/${bamName}/bamout.bam
    touch ${outputFolder}mutect2/${bamName}/f1r2.tar.gz

    $gatk --java-options "${javaOpt}" \
      GetSampleName \
        -R $refFasta \
        -I $1 \
        -O ${outputFolder}mutect2/${bamName}/${bamName}.txt \
        -encode

    $gatk --java-options "${javaOpt}" \
      Mutect2 \
        -R $refFasta \
        -I $1 \
        -tumor 'cat ${outputFolder}mutect2/${bamName}/${bamName}.txt' \
        --germline-resource $gnomad \
        -L $regions \
        -O ${outputFolder}mutect2/${bamName}/${bamName}.unfiltered.vcf \
        --bam-output ${outputFolder}mutect2/${bamName}/${bamName}.bamout.bam \
        --f1r2-tar-gz ${outputFolder}mutect2/${bamName}/f1r2.tar.gz

    $gatk --java-options "${javaOpt}" \
      GetPileupSummaries \
        -R $refFasta \
        -I $1 \
        --interval-set-rule INTERSECTION -L $regions \
        -V $variantsForContamination \
        -L $variantsForContamination \
        -O ${outputFolder}mutect2/${bamName}/${bamName}.pileups.table
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
    local inputFiles=$(printf -- "-I %s " $files)

    $gatk --java-options "${javaOpt}" \
      GatherPileupSummaries \
        --sequence-dictionary $refDict \
        $inputFiles \
        -O ${outputFolder}mutect2/gathered.pileups.table

}

# End of merging functions

# By looking into "artifact-priors.tar.gz" we can see that
# LearnReadOritentationModel defines samples correctly
# if multiple files are at the input
function learnReadOrientationModel {
    local output=$(echo ${1} | sed "s/$(basename -- ${1})/artifact-priors.tar.gz/")

    $gatk --java-options "${javaOpt}" \
      LearnReadOrientationModel \
        -I $1 \
        -O $output
}

function calculateContamination {
    local bamName=$(basename -- ${1} | cut -f 1 -d '.')

    $gatk --java-options "${javaOpt}" \
      CalculateContamination \
        -I $1 \
        -O ${outputFolder}mutect2/${bamName}/${bamName}.contamination.table \
        --tumor-segmentation ${outputFolder}mutect2/${bamName}/${bamName}.segments.table
}

function filterMutectCalls {
    local output=$(echo ${1} | sed "s/unfiltered/filtered/")
    local stats=$(echo ${1} | sed "s/unfiltered.vcf/unfiltered.vcf.stats/")
    local bamName=$(basename -- ${1} | cut -f 1 -d '.')

    $gatk --java-options "${javaOpt}" \
      FilterMutectCalls \
        -V $1 \
        -R $refFasta \
        -O $output \
        --contamination-table ${outputFolder}mutect2/${bamName}/${bamName}.contamination.table \
        --tumor-segmentation ${outputFolder}mutect2/${bamName}/${bamName}.segments.table \
        --ob-priors ${outputFolder}mutect2/${bamName}/artifact-priors.tar.gz
}

function filterAlignmentArtifacts {
    local files="${outputFolder}recalibrated/*.bam"
    # https://gatkforums.broadinstitute.org/gatk/discussion/24190/which-bam-to-give-filteralignmentartifacts

    for bam in $files; do
        local bamName=$(basename -- ${bam} | cut -f 1 -d '.')
        local vcf="${outputFolder}mutect2/${bamName}/${bamName}.filtered.vcf"
        local output=$(echo ${vcf} | sed "s/filtered/filtered.artifacts/")

        $gatk --java-options "${javaOpt}" \
          FilterAlignmentArtifacts \
            -R $refFasta \
            -V $vcf \
            -I $bam \
            --bwa-mem-index-image $refImg \
            -O $output
    done
}

makeDirectory mutect2
export -f makeDirectory
parallelRun mutect2 "${outputFolder}recalibrated/*.bam"
sleep 1

parallelRun learnReadOrientationModel "${outputFolder}mutect2/*/f1r2.tar.gz"
sleep 1

parallelRun calculateContamination "${outputFolder}mutect2/*/*.pileups.table"
sleep 1

parallelRun filterMutectCalls "${outputFolder}mutect2/*/*.unfiltered.vcf"

# There is no need to run any of these functions for all files
# in the run:
# https://github.com/broadinstitute/gatk/tree/master/scripts/mutect2_wdl
# The repository above indicates that there are no difference
# in initiation of Mutect2 pipeline functions
