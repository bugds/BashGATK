#!/bin/bash

# Multisample Mutect2 usage described here: 
# https://github.com/broadinstitute/gatk/blob/master/scripts/mutect2_wdl/mutect2_multi_sample.wdl

set -e
set -o pipefail

export outputFolder="$(realpath $1)/"
export gatk=/opt/gatk-4.1.4.1/gatk

export regions=/home/bioinfuser/NGS/Reference/intervals/2020_02_02/regions.interval_list
export refFasta=/home/bioinfuser/NGS/Reference/hg38/hg38.fasta
export refImg=/home/bioinfuser/NGS/Reference/hg38/hg38.fasta.img
export refDict=/home/bioinfuser/NGS/Reference/hg38/hg38.dict
export gnomad=/home/bioinfuser/NGS/Reference/intervals/2020_02_02/additional/AFonly.vcf
export variantsForContamination=/home/bioinfuser/NGS/Reference/intervals/2020_02_02/additional/variants_for_contamination.vcf

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
        echo $bam
        local bamName=$(basename -- ${bam} | cut -f 1 -d '.')
        makeDirectory mutect2/${bamName}

        touch ${outputFolder}mutect2/${bamName}/bamout.bam
        touch ${outputFolder}mutect2/${bamName}/f1r2.tar.gz

        $gatk --java-options "${javaOpt}" \
          GetSampleName \
            -R $refFasta \
            -I $bam \
            -O ${outputFolder}mutect2/${bamName}/${bamName}.txt \
            -encode

        $gatk --java-options "${javaOpt}" \
          Mutect2 \
            -R $refFasta \
            -I $bam \
            -tumor 'cat ${outputFolder}mutect2/${bamName}/${bamName}.txt' \
            --germline-resource $gnomad \
            -L $regions \
            -O ${outputFolder}mutect2/${bamName}/${bamName}.unfiltered.vcf \
            --bam-output ${outputFolder}mutect2/${bamName}/${bamName}.bamout.bam \
            --f1r2-tar-gz ${outputFolder}mutect2/${bamName}/f1r2.tar.gz
        
        $gatk --java-options "${javaOpt}" \
          GetPileupSummaries \
            -R $refFasta \
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
    local files="${outputFolder}mutect2/*/f1r2.tar.gz"
    local inputFiles=$(printf -- "-I %s " $files)

    $gatk --java-options "${javaOpt}" \
      LearnReadOrientationModel \
        $inputFiles \
        -O ${outputFolder}mutect2/artifact-priors.tar.gz
}

function calculateContamination {
    local files="${outputFolder}mutect2/*/*.pileups.table"

    for input in $files; do
        local bamName=$(basename -- ${input} | cut -f 1 -d '.')

        $gatk --java-options "${javaOpt}" \
          CalculateContamination \
            -I $input \
            -O ${outputFolder}mutect2/${bamName}/${bamName}.contamination.table \
            --tumor-segmentation ${outputFolder}mutect2/${bamName}/${bamName}.segments.table
    done
}

function filterMutectCalls {
    local files=${outputFolder}mutect2/*/*.unfiltered.vcf

    for vcf in $files; do
        local output=$(echo ${vcf} | sed "s/unfiltered/filtered/")
        local stats=$(echo ${vcf} | sed "s/unfiltered.vcf/unfiltered.vcf.stats/")
        local bamName=$(basename -- ${vcf} | cut -f 1 -d '.')

        $gatk --java-options "${javaOpt}" \
          FilterMutectCalls \
            -V $vcf \
            -R $refFasta \
            -O $output \
            --contamination-table ${outputFolder}mutect2/${bamName}/${bamName}.contamination.table \
            --tumor-segmentation ${outputFolder}mutect2/${bamName}/${bamName}.segments.table \
            --ob-priors ${outputFolder}mutect2/artifact-priors.tar.gz
    done
}

function filterAlignmentArtifacts {
    local files="${outputFolder}recalibrated/*.bam"
    # https://gatkforums.broadinstitute.org/gatk/discussion/24190/which-bam-to-give-filteralignmentartifacts

    for bam in $files; do
        local bamName=$(basename -- ${bam})
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

# function annotateVep {

# }

# function annotateAnnovar {

# }

# mutect2
# sleep 5
# learnReadOrientationModel
# sleep 5
# calculateContamination
# sleep 5
# filterMutectCalls
# sleep 5
fileterAlignmentArtifacts

# There is no need to run any of these functions for all files
# in the run:
# https://github.com/broadinstitute/gatk/tree/master/scripts/mutect2_wdl
# The repository above indicates that there are no difference
# in initiation of Mutect2 pipeline functions