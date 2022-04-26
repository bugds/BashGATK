#!/bin/bash

set -e
set -o pipefail

function makeDirectory {
    local newDirectory=$1

    if !( test -d "${newDirectory}/" ); then
        mkdir "${newDirectory}/"
    fi
}

function indexVcf {
    $gatk --java-options "${javaOpt}" \
      IndexFeatureFile \
        -I $gnomad
}

function selectVariants {
    makeDirectory $wd

    $gatk --java-options "${javaOpt}" \
      SelectVariants \
        -V $gnomad \
        -L $regions \
        -O ${wd}/AFonly.vcf \
        --lenient
}

function indexVcf2 {
    $gatk --java-options "${javaOpt}" \
      IndexFeatureFile \
        -I ${wd}/AFonly.vcf
}

function selectCommonBiallelicSNPs {
    $gatk --java-options "${javaOpt}" \
      SelectVariants \
        -V ${wd}/AFonly.vcf \
        -select-type SNP -restrict-alleles-to BIALLELIC \
        -select "AF > ${minimumAlleleFrequency}" \
        -O ${wd}/variants_for_contamination.vcf \
        --lenient
}

indexVcf
selectVariants
indexVcf2
selectCommonBiallelicSNPs
