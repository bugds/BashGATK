#!/bin/bash

set -e
set -o pipefail

export inputFolder=${outputFolder}/fastq
export platform=ILLUMINA
# Since we have only one lane in MiSeq
export lane=1

export regions=/home/bioinfuser/NGS/Reference/intervals/2020_02_02/capture_targets.bed
export dbSnpVcf=/home/bioinfuser/NGS/Reference/hg38/dbsnp138.vcf
export dbSnpVcfIdx=/home/bioinfuser/NGS/Reference/hg38/dbsnp138.vcf.idx
export millisVcf=/home/bioinfuser/NGS/Reference/hg38/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz
export millisVcfIdx=/home/bioinfuser/NGS/Reference/hg38/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz.tbi
export indelsVcf=/home/bioinfuser/NGS/Reference/hg38/known_indels.hg38.vcf.gz
export indelsVcfIdx=/home/bioinfuser/NGS/Reference/hg38/known_indels.hg38.vcf.gz.tbi

export bwaVersion="$($bwa 2>&1 | grep -e '^Version' | sed 's/Version: //')"
export bwaCommandline="$bwa mem -K 100000000 -p -v 3 -t 16 -Y $refFasta"
export compressionLevel=5
export parallelJobs=4

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
        echo 'No reverse fastq file for ${forwardName}'
        exit 1
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
        -O "${outputFolder}unmapped/${name}.bam" \
        -RG "rg${sample}:${lane}" \
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
    makeDirectory temporary_files

    for bam in $files; do
        gatk ValidateSamFile \
            -I $bam \
            --QUIET true \
            -M SUMMARY \
            -O "${outputFolder}temporary_files/validate.txt"
        
        local result=$(head -n 1 "${outputFolder}temporary_files/validate.txt")
        
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

    java -Dsamjdk.compression_level=${compressionLevel} ${javaOpt} -jar $picard \
      SamToFastq \
        INPUT=$1 \
        FASTQ=/dev/stdout \
        INTERLEAVE=true \
        NON_PF=true \
    | \
    ${bwaCommandline} /dev/stdin - 2> >(tee ./processing.sh_logs/${output}.bwa.stderr.log >&2) \
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
    parallel -j $parallelJobs samToFastqAndBwaMem ::: $files
}

function markDuplicates {
    local files="${outputFolder}merged/*.bam"
    local inputFiles=$(printf -- "--INPUT %s " $files)
    local javaOpt="-Xms4000m"
    makeDirectory duplicates_marked

    for bam in $files; do
        $gatk --java-options "-Dsamjdk.compression_level=${compressionLevel} ${javaOpt}" \
          MarkDuplicates \
            --INPUT $bam \
            --OUTPUT ${outputFolder}duplicates_marked/$(basename -- ${bam}) \
            --METRICS_FILE ${outputFolder}duplicates_marked/$(basename -- ${bam}).mtrx \
            --VALIDATION_STRINGENCY SILENT \
            --OPTICAL_DUPLICATE_PIXEL_DISTANCE 2500 \
            --SORTING_COLLECTION_SIZE_RATIO 0.25 \
            --ASSUME_SORT_ORDER queryname \
            --CREATE_MD5_FILE true
    done
}

function sortAndFixTags {
    local files="${outputFolder}duplicates_marked/*.bam"
    local javaOpt="-Xms4000m"
    local javaOpt2="-Xms500m"
    makeDirectory sorted

    for bam in $files; do
        base=$(basename -- "$bam")
        output=$(echo $base | cut -f 1 -d '.')

        $gatk --java-options "-Dsamjdk.compression_level=${compressionLevel} ${javaOpt}" \
          SortSam \
            --INPUT ${bam} \
            --OUTPUT /dev/stdout \
            --SORT_ORDER "coordinate" \
            --CREATE_INDEX false \
            --CREATE_MD5_FILE false \
        | \
        $gatk --java-options "-Dsamjdk.compression_level=${compressionLevel} ${javaOpt2}" \
          SetNmMdAndUqTags \
            --INPUT /dev/stdin \
            --OUTPUT ${outputFolder}sorted/${output}.bam \
            --CREATE_INDEX true \
            --CREATE_MD5_FILE true \
            --REFERENCE_SEQUENCE ${refFasta}
    done
}

function baseRecalibrator {
    local files="${outputFolder}sorted/*.bam"
    local javaOpt='-Xms4000m'
    makeDirectory temporary_files/recal_reports

    for bam in $files; do
        local bamName=$(basename -- ${bam} | cut -d "." -f 1)
        
        $gatk --java-options ${javaOpt} \
          BaseRecalibrator \
            -R $refFasta \
            -I $bam \
            --use-original-qualities \
            -O ${outputFolder}temporary_files/recal_reports/${bamName}.recal \
            --known-sites $dbSnpVcf \
            --known-sites $millisVcf \
            --known-sites $indelsVcf \
            -L $regions
    done
}

function applyBqsr {
    local files="${outputFolder}sorted/*.bam"
    makeDirectory recalibrated
    local javaOpt='-Xms3000m'
    export -f makeDirectory
    export -f applyBqsr

    for bam in $files; do
        local bamName=$(basename -- ${bam} | cut -d "." -f 1)


        $gatk --java-options ${javaOpt} \
          ApplyBQSR \
            -R $refFasta \
            -I $bam \
            # USING -L OPTION SHOULD BE AVOIDED: MATE PAIRS DISAPPEAR
            -O ${outputFolder}recalibrated/${bamName}.bam \
            -bqsr ${outputFolder}temporary_files/recal_reports/${bamName}.recal \
            --static-quantized-quals 10 --static-quantized-quals 20 --static-quantized-quals 30 \
            --add-output-sam-program-record \
            --create-output-bam-md5 \
            --use-original-qualities
    done
}

# MAIN

pairedFastQsToUnmappedBAM
sleep 5
validateSam
sleep 5
parallelMapping
sleep 5
markDuplicates
sleep 5
sortAndFixTags
sleep 5
baseRecalibrator
sleep 5
applyBqsr
