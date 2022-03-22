#!/bin/bash

set -e
set -o pipefail

function createIntervals {
    for file in ${outputFolder}sorted/*.bam
    do
        filename=$(basename -- $file)
        echo $filename
        $bedtools genomecov -bg -ibam $file \
        | awk '$4 > 100' \
        | $bedtools sort -i - \
        > ${outputFolder}${filename}_bedtools.bed
    done

    cat ${outputFolder}*_bedtools.bed \
    | $bedtools sort -i - \
    | $bedtools merge \
    > intervals.bed

    rm ${outputFolder}*_bedtools.bed
}

createIntervals
