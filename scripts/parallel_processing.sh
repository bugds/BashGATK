#!/bin/bash

set -e
set -o pipefail

function getMetadata {
    filename=$(basename -- $1)
    
    substrings=$(echo $filename | tr '_' ' ')

    for s in $substrings; do
        if [[ ${s:0:1} == $nameSubString ]]; then 
            sample=${s:1}
            name=${s:1}
        elif [[ ${s:0:1} == $laneSubString ]]; then
            lane=${s:1}
        elif [[ ${s:0:1} == "R" ]]; then
            if [[ ${s:1:2} == "1" ]]; then
                direction=true
            else
                direction=false
            fi
        fi
    done
}
export -f getMetadata

function forwardToReverse {
    local forwardPath=$1
    local forwardName=$(basename -- $1)

    local reverseName=${forwardName/_R1_/_R2_}
    local reversePath=${forwardPath/$forwardName/$reverseName}

    if test -f $reversePath; then
        reverse=$reversePath
    else
        echo 'No reverse fastq file for '${forwardName}
        exit 1
    fi
}
export -f forwardToReverse

function makeDirectory {
    local newDirectory=$1

    if !( test -d "${outputFolder}${newDirectory}/" ); then
        mkdir "${outputFolder}${newDirectory}/"
    fi
}

function fastqQualityControl {
    $fastqc -t $3 $1 -o $2
    multiqc -f $2 -o $2
}

function trimFastq {
    local files="${inputFolder}*"

    for forward in $files; do
        getMetadata $forward
        if $direction; then
            forwardToReverse $forward
            $fastp \
                -i $forward \
                -I $reverse \
                -o "${outputFolder}trimmed/$(basename -- $forward)" \
                -O "${outputFolder}trimmed/$(basename -- $reverse)" \
                -g -W 5 -q 20 -u 40 -x -3 -l 75 -c \
                -j "${outputFolder}qc_2/${output}_fastp.json" \
                -h "${outputFolder}qc_2/${output}_fastp.html" \
                --detect_adapter_for_pe \
                -w 12 $1
        fi
    done
}

function fastqToSam {
    local basename=$(basename -- $forward)

    $gatk FastqToSam \
        -F1 $forward \
        -F2 $reverse \
        -O "${outputFolder}unmapped/${basename}.bam" \
        -RG "rg${sample}:${lane}" \
        -SM $sample \
        -LB "Lib" \
        -PL $platform
}
export -f fastqToSam

function pairedFastQsToUnmappedBAM {
    forward=$1
    getMetadata $forward
    if $direction; then
        forwardToReverse $forward
        fastqToSam
    fi
}

function validateSam {
    bam=$1
    bamName=$(basename $bam)

    $gatk ValidateSamFile \
        -I $bam \
        --QUIET true \
        -M SUMMARY \
        -O "${outputFolder}temporary_files/validate${bamName}.txt"
    
    local result=$(head -n 1 "${outputFolder}temporary_files/validate${bamName}.txt")
    
    if ! [[ $result == 'No errors found' ]]; then
        echo "${bam} is invalid"
        exit 1
    fi
}

function parallelRun {  
    local files=$2

    export -f $1
    parallel -j $3 $1 ::: $files
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
        --OUTPUT "${outputFolder}premerged/${output}.bam" \
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

function mergeSamFiles {
    local javaOpt="-Xms4000m"

    filename=$(basename -- $1)
        
    substrings=$(echo $filename | tr '_' ' ')

    for s in $substrings; do
        if [[ ${s:0:1} == $nameSubString ]]; then
            name=${s:0}
        fi
    done

    if [ ! -f "${outputFolder}merged/${name}.bam" ]; then
        echo $name
        $samtools merge -r \
            ${outputFolder}merged/${name}.bam \
            ${outputFolder}premerged/*${name}*
    fi
}

function querynameSort {
    local javaOpt="-Xms4000m"
    local bamName=$(basename -- ${1} | cut -d "." -f 1)
    
    mv $1 ${1}_temp
    $samtools sort -n ${1}_temp > $1
}

function markDuplicates {
    local javaOpt="-Xms4000m"

    $gatk --java-options "-Dsamjdk.compression_level=${compressionLevel} ${javaOpt}" \
      MarkDuplicates \
        --INPUT $1 \
        --OUTPUT ${outputFolder}duplicates_marked/$(basename -- ${1}) \
        --METRICS_FILE ${outputFolder}duplicates_marked/$(basename -- ${1}).mtrx \
        --VALIDATION_STRINGENCY SILENT \
        --OPTICAL_DUPLICATE_PIXEL_DISTANCE $optimal_dup_pixel_distance \
        --SORTING_COLLECTION_SIZE_RATIO 0.25 \
        --ASSUME_SORT_ORDER queryname \
        --CREATE_MD5_FILE true
}

function sortAndFixTags {
    local javaOpt="-Xms4000m"
    local javaOpt2="-Xms500m"

    base=$(basename -- "$1")
    output=$(echo $base | cut -f 1 -d '.')

    $gatk --java-options "-Dsamjdk.compression_level=${compressionLevel} ${javaOpt}" \
      SortSam \
        --INPUT ${1} \
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
}

function baseRecalibrator {
    local javaOpt='-Xms4000m'
    local bamName=$(basename -- ${1} | cut -d "." -f 1)

    $gatk --java-options ${javaOpt} \
      BaseRecalibrator \
        -R $refFasta \
        -I $1 \
        --use-original-qualities \
        -O ${outputFolder}temporary_files/recal_reports/${bamName}.recal \
        --known-sites $dbSnpVcf \
        --known-sites $millisVcf \
        --known-sites $indelsVcf \
        -L $regions
}

function applyBqsr {
    local javaOpt='-Xms3000m'
    local bamName=$(basename -- ${1} | cut -d "." -f 1)

    $gatk --java-options ${javaOpt} \
      ApplyBQSR \
        -R $refFasta \
        -I $1 \
        -O ${outputFolder}recalibrated/${bamName}.bam \
        -bqsr ${outputFolder}temporary_files/recal_reports/${bamName}.recal \
        --static-quantized-quals 10 --static-quantized-quals 20 --static-quantized-quals 30 \
        --add-output-sam-program-record \
        --create-output-bam-md5 \
        --use-original-qualities
# USING -L OPTION SHOULD BE AVOIDED: MATE PAIRS DISAPPEAR
}

# MAIN

makeDirectory qc_1
fastqQualityControl "${inputFolder}/*" "${outputFolder}/qc_1" $parallelJobs1

makeDirectory qc_2
makeDirectory trimmed
trimFastq $parallelJobs2
sleep 1

fastqQualityControl "${outputFolder}trimmed/*" "${outputFolder}/qc_2" $parallelJobs1
sleep 1

makeDirectory unmapped
parallelRun pairedFastQsToUnmappedBAM "${outputFolder}trimmed/*" $parallelJobs3
sleep 1
rm -r "${outputFolder}trimmed"

makeDirectory temporary_files
parallelRun validateSam "${outputFolder}unmapped/*" $parallelJobs4
sleep 1

makeDirectory unmerged
makeDirectory premerged
parallelRun samToFastqAndBwaMem "${outputFolder}unmapped/*" $parallelJobs5
sleep 1
rm -r "${outputFolder}unmapped"

makeDirectory merged
parallelRun mergeSamFiles "${outputFolder}premerged/*" $parallelJobs6
sleep 1
rm -r "${outputFolder}unmerged"
rm -r "${outputFolder}premerged"

parallelRun querynameSort "${outputFolder}merged/*.bam" $parallelJobs7
sleep 1
makeDirectory duplicates_marked
parallelRun markDuplicates "${outputFolder}merged/*.bam" $parallelJobs8
sleep 1

makeDirectory sorted
parallelRun sortAndFixTags "${outputFolder}duplicates_marked/*.bam" $parallelJobs9
sleep 1

makeDirectory temporary_files/recal_reports
parallelRun baseRecalibrator "${outputFolder}sorted/*.bam" $parallelJobs10
sleep 1

makeDirectory recalibrated
parallelRun applyBqsr "${outputFolder}sorted/*.bam" $parallelJobs11
