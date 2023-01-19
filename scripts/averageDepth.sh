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
    
    echo "\nTotal reads:" >> ${outputFolder}depths/${bamName}.txt
    samtools view -c $1 >> ${outputFolder}depths/${bamName}.txt
    echo "\nMapped reads:" >> ${outputFolder}depths/${bamName}.txt
    $samtools view -c -F 260 $1 >> ${outputFolder}depths/${bamName}.txt
}

makeDirectory depths
parallelRun getDepths "${outputFolder}sorted/*.bam"
cat ${outputFolder}depths/*.txt > "${outputFolder}depths/depthReport.out"
