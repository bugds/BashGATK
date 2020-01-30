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
export refDict=/home/bioinfuser/NGS/Reference/hg38/hg38.dict
export dbSnpVcf=/home/bioinfuser/NGS/Reference/hg38/dbsnp138.vcf
export dbSnpVcfIdx=/home/bioinfuser/NGS/Reference/hg38/dbsnp138.vcf.idx
export millisVcf=/home/bioinfuser/NGS/Reference/hg38/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz
export millisVcfIdx=/home/bioinfuser/NGS/Reference/hg38/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz.tbi
export indelsVcf=/home/bioinfuser/NGS/Reference/hg38/known_indels.hg38.vcf.gz
export indelsVcfIdx=/home/bioinfuser/NGS/Reference/hg38/known_indels.hg38.vcf.gz.tbi

export bwaVersion="$($bwa 2>&1 | grep -e '^Version' | sed 's/Version: //')"
export bwaCommandline="$bwa mem -K 100000000 -p -v 3 -t 16 -Y $refFasta"
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
        -O "${outputFolder}unmapped/${name}.bam" \
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
    ${bwaCommandline} /dev/stdin - 2> >(tee ./somatic.sh_logs/${output}.bwa.stderr.log >&2) \
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
    local files="${outputFolder}merged/*.bam"
    # local inputFiles=$(printf -- "--INPUT %s " $files)
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
            --ASSUME_SORT_ORDER "queryname" \
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

function createSequenceGroupingTSV {
    python - << EOF
with open("${refDict}", "r") as ref_dict_file:
    sequence_tuple_list = []
    longest_sequence = 0
    for line in ref_dict_file:
        if line.startswith("@SQ"):
            line_split = line.split("\t")
            # (Sequence_Name, Sequence_Length)
            sequence_tuple_list.append((line_split[1].split("SN:")[1], int(line_split[2].split("LN:")[1])))
    longest_sequence = sorted(sequence_tuple_list, key=lambda x: x[1], reverse=True)[0][1]
# We are adding this to the intervals because hg38 has contigs named with embedded colons (:) and a bug in 
# some versions of GATK strips off the last element after a colon, so we add this as a sacrificial element.
hg38_protection_tag = ":1+"
# initialize the tsv string with the first sequence
tsv_string = sequence_tuple_list[0][0] + hg38_protection_tag
temp_size = sequence_tuple_list[0][1]
for sequence_tuple in sequence_tuple_list[1:]:
    if temp_size + sequence_tuple[1] <= longest_sequence:
        temp_size += sequence_tuple[1]
        tsv_string += "\t" + sequence_tuple[0] + hg38_protection_tag
    else:
        tsv_string += "\n" + sequence_tuple[0] + hg38_protection_tag
        temp_size = sequence_tuple[1]
# add the unmapped sequences as a separate line to ensure that they are recalibrated as well
with open("${outputFolder}temporary_files/sequence_grouping.txt","w") \
  as tsv_file:
    tsv_file.write(tsv_string)
    tsv_file.close()

tsv_string += '\n' + "unmapped"

with open("${outputFolder}temporary_files/sequence_grouping_with_unmapped.txt","w") \
  as tsv_file_with_unmapped:
    tsv_file_with_unmapped.write(tsv_string)
    tsv_file_with_unmapped.close()
EOF
}

function baseRecalibrator {
    local group=$(cat /dev/stdin)
    local intervals=$(printf -- "-L %s " $group)
    local knownSitesInput=$(printf -- "--known-sites %s " $knownSites)
    local javaOpt='-Xms4000m'
    local bamName=$(basename -- ${bam} | cut -d "." -f 1)
    makeDirectory temporary_files/${bamName}_recalibration_report
    local mark=$(echo $group[0] | cut -d ' ' -f 1)

    $gatk --java-options ${javaOpt} \
      BaseRecalibrator \
        -R $refFasta \
        -I $bam \
        --use-original-qualities \
        -O ${outputFolder}temporary_files/${bamName}_recalibration_report/${mark}.txt \
        --known-sites $dbSnpVcf \
        --known-sites $millisVcf \
        --known-sites $indelsVcf \
        $intervals
}

function parallelRecalibration {
    # Reading to an array:
    #readarray -t seqGroup <${outputFolder}temporary_files/sequence_grouping.txt
    # Take all elements of an array
    #parallel baseRecalibrator ::: "${seqGroup[@]}"
    local files="${outputFolder}sorted/*.bam"

    export -f makeDirectory
    export -f baseRecalibrator
    for bam in $files; do
        export bam=$bam
        # Had to use a pipe in GNU parallel due to problems with passing large
        # strings to parallel directly (N1 means 1 line at a time)
        cat ${outputFolder}temporary_files/sequence_grouping.txt | parallel -N1 --pipe baseRecalibrator
    done
}

# MAIN

# AlreadyDone:
# pairedFastQsToUnmappedBAM
# validateSam
# parallelMapping
# sortAndFixTags
# markDuplicates
# sortAndFixTags
# createSequenceGroupingTSV - needs to be done each reference update

parallelRecalibration

# To do:
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
