#!/bin/bash

inputFolder="$(realpath $1)/"
outputFolder="$(realpath $2)/"
platform=ILLUMINA
runNumber=$3
gatk=/opt/gatk-4.1.4.1/gatk
samtools=/opt/gatk4-data-processing/samtools-1.3.1/samtools
picard=/opt/gatk4-data-processing/picard-2.16.0/picard.jar
bwa=/opt/gatk4-data-processing/bwa-0.7.15/bwa
refFasta=/home/bioinfuser/NGS/Reference/hg38/hg38.fasta

bwaVersion="$($bwa 2>&1 | grep -e '^Version' | sed 's/Version: //')"
bwaCommandline=" mem -K 100000000 -p -v 3 -t 16 -Y $refFasta"

# USE THESE VARIABLES WITH $ ONLY

function getMetadata {
    filename=$(basename -- $1)
    
    substrings=$(echo $filename | tr '_' ' ')

    for s in $substrings; do
        if [[ ${s:0:1} == "S" ]]; then 
            sample=${s:1}
        elif [[ ${s:0:3} == "MED" ]]; then
            name=$s
        elif [[ ${s:0:1} == "L" ]]; then
            library=${s:1}
        elif [[ ${s:0:1} == "R" ]]; then
            if [[ ${s:1:2} == "1" ]]; then
                direction=true
            else
                direction=false
            fi
        fi
    done
}

function forwardToReverse {
    local forwardPath=$1
    local forwardName=$(basename -- $1)

    local reverseName=${forwardName/_R1_/_R2_}
    local reversePath=${forwardPath/$forwardName/$reverseName}

    if test -f $reversePath; then
        reverse=$reversePath
    else
        echo 'ERROR'
    fi
}

function makeDirectory {
    local newDirectory=$1

    if !( test -d "${outputFolder}${newDirectory}/" ); then
        mkdir "${outputFolder}${newDirectory}/"
    fi
}

function fastqToSam {
    gatk FastqToSam \
        -F1 $forward \
        -F2 $reverse \
        -O "${outputFolder}unmapped/${name}.unmapped.bam" \
        -RG $runNumber \
        -SM $sample \
        -LB $library \
        -PL $platform
}

function pairedFastQsToUnmappedBAM {
    local files="${inputFolder}*"
    makeDirectory unmapped

    for forward in $files; do
        getMetadata $forward
        if $direction; then
            forwardToReverse $forward
            fastqToSam
        fi
    done
}

function validateSam {
    local files="${outputFolder}unmapped/*"
    makeDirectory temporaryFiles

    for bam in $files; do
        gatk ValidateSamFile \
            -I $bam \
            --QUIET true \
            -M SUMMARY \
            -O "${outputFolder}temporaryFiles/validate.txt"
        
        local result=$(head -n 1 "${outputFolder}temporaryFiles/validate.txt")
        
        if ! [[ $result == 'No errors found' ]]; then
            echo "${bam} is invalid"
            exit 1
        fi
    done
}

function samToFastqAndBwaMem {
    local bam=$(basename -- "$1")
    local output=$(echo "$bam" | cut -f 1 -d '.')
    
    java -jar $picard SamToFastq \
        INPUT=$1 \
        FASTQ=/dev/stdout \
        INTERLEAVE=true \
        NON_PF=true \
    | \
    ${bwa}${bwaCommandline} /dev/stdin -  2> >(tee ./${0}_logs/${output}.bwa.stderr.log >&2) \
    | \
    $samtools view -1 - > ${outputFolder}unmerged/${output}.bam
}

function parallelMapping {
      
    local files="${outputFolder}unmapped/*"
    makeDirectory unmerged
    export -f samToFastqAndBwaMem
    parallel samToFastqAndBwaMem ::: $files
}

# MAIN

# AlreadyDone:
# pairedFastQsToUnmappedBAM
# validateSam

parallelMapping

# To do:
# scatter:
#   samToFastqAndBwaMem: picard samToFastq, bwa, samtools
#   MergeBamAlignment
# MarkDuplicates
# SortAndFixTags
# CreateSequenceGroupingTSV
# scatter:
#   BaseRecalibrator
# GatherBqsrReports
# scatter:
#   ApplyBQSR
# GatherBamFiles


# ValidateSamFile

# gatk4-data-processing
# SamToFastqAndBwaMem
# MergeBamAlignment
# MarkDuplicates
# SortAndFixTags
# CreateSequenceGroupingTSV
#
