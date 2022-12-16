#!/bin/bash

set -e
set -o pipefail

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

function getMetadata {
    filename=$(basename -- $1)
    
    substrings=$(echo $filename | tr '_' ' ')

    for s in $substrings; do
        if [[ ${s:0:1} == $nameSubString ]]; then 
            sample=${s:1}
            name=${s:1}
        # elif [[ ${s:0:3} == "MED" ]]; then
        #    name=$s
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

function pairedFastQsToUnmappedBAM {
    local files="${outputFolder}trimmed/*"

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

    for bam in $files; do
        $gatk ValidateSamFile \
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

function extractUMIs {
    local bam=$(basename -- "$1")
    local output=$(echo $bam | cut -f 1 -d '.')
    local javaOpt="-Xms3000m"

    java -Dsamjdk.compression_level=${compressionLevel} ${javaOpt} -jar $fgbio ExtractUmisFromBam \
        -i $1 \
        -r 3M3S+T 3M3S+T \
        -t RX \
        -a true \
        -o "${outputFolder}UMIs_extracted/${output}.ubam"
}

function samToFastq {
    $gatk SamToFastq \
        -i $1 \
        -F $2 \
        -F2 $3 \
        --CLIPPING_ATTRIBUTE XT \
        --CLIPPING_ACTION 2
}

function qcAndAlign {
    local bam=$(basename -- "$1")
    local output=$(echo $bam | cut -f 1 -d '.')
    local javaOpt="-Xms3000m"

    samToFastq \
        $1 \
        "${outputFolder}fastq_for_qc/${output}_R1.fastq" \
        "${outputFolder}fastq_for_qc/${output}_R2.fastq"
    
    $fastp \
        -i "${outputFolder}fastq_for_qc/${output}_R1.fastq" \
        -o "${outputFolder}fastq_trimmed/${output}_R1.fastq" \
        -I "${outputFolder}fastq_for_qc/${output}_R2.fastq" \
        -O "${outputFolder}fastq_trimmed/${output}_R2.fastq" \
        -g -W 5 -q 20 -u 40 -x -3 -l 75 -c \
        -j "${outputFolder}fastq_trimmed/${output}_fastp.json" \
        -h "${outputFolder}fastq_trimmed/${output}_fastp.html" \
        -w 12

    ${bwaCommandline1} \
        "${outputFolder}fastq_trimmed/${output}_R1.fastq" \
        "${outputFolder}fastq_trimmed/${output}_R2.fastq" \
    | \
    $samtools view -Sb > "${outputFolder}unmerged/${output}.bam"

    $gatk --java-options "-Dsamjdk.compression_level=${compressionLevel} ${javaOpt}" \
      MergeBamAlignment \
        --ATTRIBUTES_TO_RETAIN X0 \
        --ATTRIBUTES_TO_REMOVE NM \
        --ATTRIBUTES_TO_REMOVE MD \
        --ALIGNED_BAM "${outputFolder}unmerged/${output}.bam" \
        --UNMAPPED_BAM ${1} \
        --OUTPUT "${outputFolder}premerged/${output}.bam" \
        --REFERENCE_SEQUENCE ${refFasta} \
        --SORT_ORDER "queryname" \
        --ALIGNED_READS_ONLY true \
        --MAX_INSERTIONS_OR_DELETIONS -1 \
        --PRIMARY_ALIGNMENT_STRATEGY MostDistant \
        --ALIGNER_PROPER_PAIR_FLAGS true \
        --CLIP_OVERLAPPING_READS false
}

function mergeSamFiles {
    local javaOpt="-Xms4000m"
    local files="${outputFolder}premerged/*"

    for file in $files; do
        filename=$(basename -- $file)
        
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
    done
}

function onlyProperPairs {
    local bam=$(basename -- "$1")
    local output=$(echo $bam | cut -f 1 -d '.')

    $samtools view -f 2 \
        -bh $1\
        > "${outputFolder}only_proper/${output}.bam"
}

function groupSameUMIReads {
    local bam=$(basename -- "$1")
    local output=$(echo $bam | cut -f 1 -d '.')
    local javaOpt="-Xms3000m"

    java -Dsamjdk.compression_level=${compressionLevel} ${javaOpt} -jar $fgbio GroupReadsByUmi \
        --input=$1 \
        --output="${outputFolder}reads_grouped/${output}.bam" \
        --strategy=adjacency \
        --edits=1 \
        -t RX \
        -f "${outputFolder}reads_grouped/${output}.txt"

    java -Dsamjdk.compression_level=${compressionLevel} ${javaOpt} -jar $fgbio CallMolecularConsensusReads \
        --input="${outputFolder}reads_grouped/${output}.bam" \
        --output="${outputFolder}reads_called/${output}.bam" \
        --error-rate-post-umi 40 \
        --error-rate-pre-umi 45 \
        --output-per-base-tags false \
        --min-reads 2 \
        --max-reads 50 \
        --min-input-base-quality 20 \
        --read-name-prefix="consensus"
}

function realign {
    local bam=$(basename -- "$1")
    local output=$(echo $bam | cut -f 1 -d '.')
    local javaOpt="-Xms3000m"

    samToFastq \
        $1 \
        "${outputFolder}fastq_for_realignment/${output}_R1.fastq" \
        "${outputFolder}fastq_for_realignment/${output}_R2.fastq"

    ${bwaCommandline2} \
        "${outputFolder}fastq_for_realignment/${output}_R1.fastq" \
        "${outputFolder}fastq_for_realignment/${output}_R2.fastq" \
    | \
    $samtools view -bh > "${outputFolder}realigned/${output}.bam"
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

makeDirectory unmapped
pairedFastQsToUnmappedBAM
sleep 1

makeDirectory temporary_files
validateSam
sleep 1

makeDirectory UMIs_extracted
parallelRun extractUMIs "${outputFolder}unmapped/*"

makeDirectory fastq_for_qc
makeDirectory fastq_trimmed
makeDirectory unmerged
makeDirectory premerged
parallelRun qcAndAlign "${outputFolder}UMIs_extracted/*"
sleep 1

makeDirectory merged
mergeSamFiles
sleep 1

makeDirectory only_proper
parallelRun onlyProperPairs "${outputFolder}merged/*"
sleep 1

makeDirectory reads_grouped
makeDirectory reads_called
parallelRun groupSameUMIReads "${outputFolder}only_proper/*"

makeDirectory fastq_for_realignment
makeDirectory realigned
parallelRun realign "${outputFolder}reads_called/*"

makeDirectory sorted
parallelRun sortAndFixTags "${outputFolder}realigned/*.bam"
sleep 1

makeDirectory temporary_files/recal_reports
parallelRun baseRecalibrator "${outputFolder}sorted/*.bam"
sleep 1

makeDirectory recalibrated
parallelRun applyBqsr "${outputFolder}sorted/*.bam"
