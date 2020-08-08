#!/bin/bash

set -e
set -o pipefail

export cnvkit='python3 /home/bioinfuser/NGS/Software/cnvkit/cnvkit.py'
export regions='/home/bioinfuser/NGS/Reference/intervals/2020_02_02/cnvkit_targets.bed'
export antitarget='/home/bioinfuser/NGS/Reference/intervals/2020_02_02/cnvkit_antitargets.bed'
export mappable='/home/bioinfuser/NGS/Reference/hg38/access-hg38.bed'
export refFlat='/home/bioinfuser/NGS/Reference/hg38/hg38.refFlat.txt'
export AFonly='/home/bioinfuser/NGS/Reference/intervals/2020_02_02/additional/AFonly.vcf'
export LC_NUMERIC='en_US.UTF-8'

export javaOpt="-Xms3000m"

function makeDirectory {
    local newDirectory=$1

    if !( test -d "${outputFolder}/${newDirectory}/" ); then
        mkdir "${outputFolder}/${newDirectory}/"
    fi
}

function putTimestamp {
    # To ensure timestamp on bai later than on bam;
    # or else cnvkit will try to modify bai!!!

    for f in $outputFolder/recalibrated/*.bai; do
        touch $f
    done
    for f in $outputFolder/normals/*.bai; do
        touch $f
    done
}

function selectGermlineVcf {
    for vcf in $outputFolder/deepvariant/*.vcf.gz; do
        baseName=$(basename -- $vcf | cut -d '.' -f1)

        gatk --java-options "${javaOpt}" \
          SelectVariants \
            -V $vcf \
            --concordance $AFonly \
            -O $outputFolder/deepvariant/$baseName.baf.vcf
    done
}

function cnvKitAutobin {
    #$cnvkit access \
    #    $refFasta \
    #    -o $mappable
    
    cd $outputFolder/cnvkit/
    # IDK where autobin generates files, probably in the wd
    
    $cnvkit autobin \
        $outputFolder/{normals,recalibrated}/*.bam \
        -t $regions \
        -g $mappable \
        --annotate $refFlat \
        --short-names
}

function cnvKitCoverage {
    for bam in $outputFolder/normals/*.bam; do
        base=$(basename -- $bam | cut -d . -f 1)
        
        $cnvkit coverage \
            $bam \
            $outputFolder/cnvkit/cnvkit_targets.target.bed \
            -o $outputFolder/cnvkit/normal/$base.targetcoverage.cnn
            
        $cnvkit coverage \
            $bam \
            $outputFolder/cnvkit/cnvkit_targets.antitarget.bed \
            -o $outputFolder/cnvkit/normal/$base.antitargetcoverage.cnn
    done

    for bam in $outputFolder/recalibrated/*.bam; do
        base=$(basename -- $bam | cut -d . -f 1)
        
        $cnvkit coverage \
            $bam \
            $outputFolder/cnvkit/cnvkit_targets.target.bed \
            -o $outputFolder/cnvkit/pathologic/$base.targetcoverage.cnn
            
        $cnvkit coverage \
            $bam \
            $outputFolder/cnvkit/cnvkit_targets.antitarget.bed \
            -o $outputFolder/cnvkit/pathologic/$base.antitargetcoverage.cnn
    done
}

function cnvKitReference {
    $cnvkit reference \
        $outputFolder/cnvkit/normal/*.{,anti}targetcoverage.cnn \
        --fasta $refFasta \
        -o $outputFolder/cnvkit/normal/reference.cnn
}

function cnvKitCall {
    # -m hybrid: SeqCap hybridization enrichment was performed
    # https://github.com/bcbio/bcbio-nextgen/issues/2171 - drop-low-coverage
    
    for bam in $outputFolder/recalibrated/*.bam; do
        base=$(basename -- $bam | cut -d . -f 1)
        
        $cnvkit fix \
            $outputFolder/cnvkit/pathologic/$base.targetcoverage.cnn \
            $outputFolder/cnvkit/pathologic/$base.antitargetcoverage.cnn \
            $outputFolder/cnvkit/normal/reference.cnn \
            -o $outputFolder/cnvkit/pathologic/$base.cnr
            
        $cnvkit segment \
            $outputFolder/cnvkit/pathologic/$base.cnr \
            --drop-low-coverage \
            -t 1e-6 \
            -o $outputFolder/cnvkit/pathologic/$base.cns
        
        cd $outputFolder/cnvkit/pathologic
            
        $cnvkit segmetrics \
            $outputFolder/cnvkit/pathologic/$base.cnr \
            -s $outputFolder/cnvkit/pathologic/$base.cns \
            --ci
            
        $cnvkit call \
            $outputFolder/cnvkit/pathologic/$base.segmetrics.cns \
            --center median \
            --filter ci \
            -v $outputFolder/deepvariant/$base.baf.vcf \
            -o $outputFolder/cnvkit/pathologic/$base.vcf.cns
            #--purity $(cat $outputFolder/purity/$base.txt) \
            #--purity $(cat $outputFolder/doublePurity/$base.txt) \
            
        $cnvkit scatter \
            $outputFolder/cnvkit/pathologic/$base.cnr \
            -s $outputFolder/cnvkit/pathologic/$base.vcf.cns \
            -o $outputFolder/cnvkit/pathologic/$base-scatter.png
    done
}


makeDirectory cnvkit
makeDirectory cnvkit/normal
makeDirectory cnvkit/pathologic

#putTimestamp
#selectGermlineVcf

cnvKitAutobin
cnvKitCoverage
cnvKitReference
cnvKitCall
