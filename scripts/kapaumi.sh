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
    presubstrings=$(echo $filename | cut -d '.' -f1)
    substrings=$(echo $presubstrings | tr '_' ' ')

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
        -o "${outputFolder}UMIs_extracted/${output}.bam"
}

function qcAndAlign {
    local bam=$(basename -- "$1")
    local output=$(echo $bam | cut -f 1 -d '.')
    local javaOpt="-Xms3000m"

    inputSam=$1

    $gatk SamToFastq \
        -I $inputSam \
        -F "${outputFolder}fastq_for_qc/${output}_R1.fastq.gz" \
        -F2 "${outputFolder}fastq_for_qc/${output}_R2.fastq.gz" \
        --CLIPPING_ATTRIBUTE XT \
        --CLIPPING_ACTION 2
    
    $fastp \
        -i "${outputFolder}fastq_for_qc/${output}_R1.fastq.gz" \
        -o "${outputFolder}fastq_trimmed/${output}_R1.fastq.gz" \
        -I "${outputFolder}fastq_for_qc/${output}_R2.fastq.gz" \
        -O "${outputFolder}fastq_trimmed/${output}_R2.fastq.gz" \
        -g -W 5 -q 20 -u 40 -x -3 -l 75 -c \
        -j "${outputFolder}fastq_trimmed/${output}_fastp.json" \
        -h "${outputFolder}fastq_trimmed/${output}_fastp.html" \
        --adapter_sequence=AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
        --adapter_sequence_r2=AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
        -w 12
    
    rm "${outputFolder}fastq_for_qc/${output}_R1.fastq.gz"
    rm "${outputFolder}fastq_for_qc/${output}_R2.fastq.gz"

    ${bwaCommandline1} \
        "${outputFolder}fastq_trimmed/${output}_R1.fastq.gz" \
        "${outputFolder}fastq_trimmed/${output}_R2.fastq.gz" \
    | \
    $samtools view -Sb > "${outputFolder}unmerged/${output}.bam"

    rm "${outputFolder}fastq_trimmed/${output}_R1.fastq.gz"
    rm "${outputFolder}fastq_trimmed/${output}_R2.fastq.gz"

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
  
    inputSam=$1

    $gatk SamToFastq \
        -I $inputSam \
        -F "${outputFolder}fastq_for_realignment/${output}_R1.fastq" \
        -F2 "${outputFolder}fastq_for_realignment/${output}_R2.fastq" \
        --CLIPPING_ATTRIBUTE XT \
        --CLIPPING_ACTION 2

    ${bwaCommandline2} \
        "${outputFolder}fastq_for_realignment/${output}_R1.fastq" \
        "${outputFolder}fastq_for_realignment/${output}_R2.fastq" \
    | \
    $samtools view -bh > "${outputFolder}realigned/${output}.bam"
    
    $gatk --java-options "-Dsamjdk.compression_level=${compressionLevel} ${javaOpt}" \
      SortSam \
        --INPUT "${outputFolder}reads_called/${output}.bam" \
        --OUTPUT "${outputFolder}reads_called/${output}_sorted.bam" \
        --SORT_ORDER "queryname" \
        --CREATE_INDEX false \
        --CREATE_MD5_FILE false \

    $gatk --java-options "-Dsamjdk.compression_level=${compressionLevel} ${javaOpt}" \
      MergeBamAlignment \
        --ATTRIBUTES_TO_RETAIN X0 \
        --ATTRIBUTES_TO_REMOVE NM \
        --ATTRIBUTES_TO_REMOVE MD \
        --ALIGNED_BAM "${outputFolder}realigned/${output}.bam" \
        --UNMAPPED_BAM "${outputFolder}reads_called/${output}_sorted.bam" \
        --OUTPUT "${outputFolder}realigned_merged/${output}.bam" \
        --REFERENCE_SEQUENCE ${refFasta} \
        --SORT_ORDER "queryname" \
        --ALIGNED_READS_ONLY true \
        --MAX_INSERTIONS_OR_DELETIONS -1 \
        --PRIMARY_ALIGNMENT_STRATEGY MostDistant \
        --ALIGNER_PROPER_PAIR_FLAGS true \
        --CLIP_OVERLAPPING_READS false
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
        --OPTICAL_DUPLICATE_PIXEL_DISTANCE 2500 \
        --SORTING_COLLECTION_SIZE_RATIO 0.25 \
        --ASSUME_SORT_ORDER queryname \
        --CREATE_MD5_FILE true
}

function getDepths {
    local bamName=$(basename -- ${1} | cut -f 1 -d '.')
    echo $bamName > ${outputFolder}depths/${bamName}.txt

    $samtools depth \
        -b $regions \
        $1 \
        -a \
        | \
        awk '$3<10{$2=1}1' \
        | \
        awk '$3>=10{$2=0}1' \
        | \
        awk '{sum+=$3; sumsq+=$3*$3; subsum+=$2} END {print "Average = ",sum/NR; print "Stdev = ",sqrt(sumsq/NR - (sum/NR)**2); print "Sub-ten = ",subsum/NR; print "Total nucleotides: ",NR}' \
        >> ${outputFolder}depths/${bamName}.txt
    
    echo "Total reads:" >> ${outputFolder}depths/${bamName}.txt
    $samtools view -c $1 >> ${outputFolder}depths/${bamName}.txt
    echo "Mapped reads:" >> ${outputFolder}depths/${bamName}.txt
    $samtools view -c -F 260 $1 >> ${outputFolder}depths/${bamName}.txt
}

# MAIN

makeDirectory unmapped
parallelRun pairedFastQsToUnmappedBAM "${outputFolder}fastq/*"
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

rm -r "${outputFolder}fastq_for_qc"
rm -r "${outputFolder}fastq_trimmed"
rm -r "${outputFolder}unmerged"
rm -r "${outputFolder}UMIs_extracted"
sleep 1

makeDirectory merged
parallelRun mergeSamFiles "${outputFolder}premerged/*"
sleep 1

rm -r "${outputFolder}unmapped"
rm -r "${outputFolder}premerged"
sleep 1

### SOMATIC PART

makeDirectory only_proper
parallelRun onlyProperPairs "${outputFolder}merged/*"
sleep 1

makeDirectory reads_grouped
makeDirectory reads_called
parallelRun groupSameUMIReads "${outputFolder}only_proper/*"

rm -r "${outputFolder}only_proper"
sleep 1

makeDirectory fastq_for_realignment
makeDirectory realigned
makeDirectory realigned_merged
parallelRun realign "${outputFolder}reads_called/*"

rm -r "${outputFolder}reads_grouped"
rm -r "${outputFolder}reads_called"
sleep 1

makeDirectory sorted
parallelRun sortAndFixTags "${outputFolder}realigned_merged/*.bam"
sleep 1

rm -r "${outputFolder}fastq_for_realignment"
rm -r "${outputFolder}realigned"
rm -r "${outputFolder}realigned_merged"
sleep 1

makeDirectory temporary_files/recal_reports
parallelRun baseRecalibrator "${outputFolder}sorted/*.bam"
sleep 1

makeDirectory recalibrated
parallelRun applyBqsr "${outputFolder}sorted/*.bam"

makeDirectory depths
parallelRun getDepths "${outputFolder}sorted/*.bam"
cat ${outputFolder}depths/*.txt > "${outputFolder}depths/depthReport_somatic.out"

rm -r "${outputFolder}sorted"
sleep 1

### GERMINAL PART

parallelRun querynameSort "${outputFolder}merged/*.bam"
sleep 1
makeDirectory duplicates_marked
parallelRun markDuplicates "${outputFolder}merged/*.bam"
sleep 1

rm -r "${outputFolder}merged"
sleep 1

makeDirectory sorted
parallelRun sortAndFixTags "${outputFolder}duplicates_marked/*.bam"
sleep 1

rm -r "${outputFolder}duplicates_marked"
sleep 1

parallelRun getDepths "${outputFolder}sorted/*.bam"
cat ${outputFolder}depths/*.txt > "${outputFolder}depths/depthReport_germinal.out"

sleep 1