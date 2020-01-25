#!/bin/bash

export inputFolder="$(realpath $1)/"
export outputFolder="$(realpath $2)/"
export platform=ILLUMINA
export runNumber=$3
export gatk=/opt/gatk-4.1.4.1/gatk
export samtools=/opt/gatk4-data-processing/samtools-1.3.1/samtools
export picard=/opt/gatk4-data-processing/picard-2.16.0/picard.jar
export bwa=/opt/gatk4-data-processing/bwa-0.7.15/bwa
export refFasta=/home/bioinfuser/NGS/Reference/hg38/hg38.fasta

export bwaVersion="$($bwa 2>&1 | grep -e '^Version' | sed 's/Version: //')"
export bwaCommandline=" mem -K 100000000 -p -v 3 -t 16 -Y $refFasta"
export compressionLevel=5

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
    local output=$(echo $bam | cut -f 1 -d '.')
    local javaOpt="-Xms3000m"

    java --java-options "-Dsamjdk.compression_level=${compressionLevel} ${javaOpt}" -jar $picard \
      SamToFastq \
        INPUT=$1 \
        FASTQ=/dev/stdout \
        INTERLEAVE=true \
        NON_PF=true \
    | \
    ${bwa}${bwaCommandline} /dev/stdin - 2> >(tee ./somatic.sh_logs/${output}.bwa.stderr.log >&2) \
    | \
    $samtools view -1 - > ${outputFolder}unmerged/${output}.bam
    
    $gatk --java-options "-Dsamjdk.compression_level=${compressionLevel} ${javaOpt}" \
      MergeBamAlignment \
        --VALIDATION_STRINGENCY SILENT \
        --EXPECTED_ORIENTATIONS FR \
        --ATTRIBUTES_TO_RETAIN X0 \
        --ALIGNED_BAM "${outputFolder}unmerged/${output}.bam" \
        --UNMAPPED_BAM ${1} \
        --OUTPUT "${outputFolder}merged/${output}.bam" \
        --REFERENCE_SEQUENCE ${refFasta} \
        --PAIRED_RUN true \
        --SORT_ORDER "unsorted" \
        --IS_BISULFITE_SEQUENCE false \
        --ALIGNED_READS_ONLY false \
        --CLIP_ADAPTERS false \
        --MAX_RECORDS_IN_RAM 2000000 \
        --ADD_MATE_CIGAR true \
        --MAX_INSERTIONS_OR_DELETIONS -1 \
        --PRIMARY_ALIGNMENT_STRATEGY MostDistant \
        --PROGRAM_RECORD_ID "bwamem" \
        --PROGRAM_GROUP_VERSION "${bwaVersion}" \
        --PROGRAM_GROUP_COMMAND_LINE "${bwaCommandline}" \
        --PROGRAM_GROUP_NAME "bwamem" \
        --UNMAPPED_READ_STRATEGY COPY_TO_TAG \
        --ALIGNER_PROPER_PAIR_FLAGS true \
        --UNMAP_CONTAMINANT_READS true
}

function parallelMapping {  
    local files="${outputFolder}unmapped/*"

    makeDirectory unmerged
    makeDirectory merged
    export -f samToFastqAndBwaMem
    parallel samToFastqAndBwaMem ::: $files
}

function markDuplicates {
    local files="${outputFolder}merged/*"
    local inputFiles=$(printf -- "--INPUT %s " $files)
    local javaOpt="-Xms4000m"
    makeDirectory duplicates_marked

    $gatk --java-options "-Dsamjdk.compression_level=${compressionLevel} ${javaOpt}" \
      MarkDuplicates \
        --INPUT ${inputFiles}\
        --OUTPUT ${outputFolder}duplicates_marked/$(basename -- ${outputFolder}).marked.bam \
        --METRICS_FILE ${outputFolder}duplicates_marked/$(basename -- ${outputFolder}).mtrx \
        --VALIDATION_STRINGENCY SILENT \
        --OPTICAL_DUPLICATE_PIXEL_DISTANCE 2500 \
        --ASSUME_SORT_ORDER "queryname" \
        --CREATE_MD5_FILE true
}

# MAIN

# AlreadyDone:
# pairedFastQsToUnmappedBAM
# validateSam
# parallelMapping

markDuplicates

# To do:
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
