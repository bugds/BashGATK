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
        awk '{sum+=$3; sumsq+=$3*$3} END { print "Average = ",sum/NR; print "Stdev = ",sqrt(sumsq/NR - (sum/NR)**2)}' \
        >> ${outputFolder}depths/${bamName}.txt
}

makeDirectory depths
parallelRun getDepths "${outputFolder}sorted/*.bam"
cat ${outputFolder}depths/*.txt > "${outputFolder}depths/depthReport.out"
