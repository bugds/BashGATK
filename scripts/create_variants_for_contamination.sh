#!/bin/bash

set -e
set -o pipefail

# Given a gnomAD vcf, produce:
# 1) a sites-only vcf where the only INFO field is allele frequency (AF)
#   this is used as the gnomad input to mutect2.wdl
# 2) this sites-only vcf further restricted to a given minimum allele frequency
#   this is used as the variants_for_contamination input to mutect2.wdl

export wd=/home/bioinfuser/NGS/Reference/intervals/2020_02_02/additional/
# export outputFolder=/home/bioinfuser/NGS/Reference/hg38/Mutect2/
export regions=/home/bioinfuser/NGS/Reference/intervals/2020_02_02/regions.bed

export gnomad=/home/bioinfuser/NGS/Reference/hg38/Mutect2/gnomad.exomes.r2.1.1.sites.liftover_grch38.vcf
export variantsForContamination=/home/bioinfuser/NGS/Reference/hg38/Mutect2/variants_for_contamination.vcf
export minimumAlleleFrequency=0.05

export javaOpt="-Xms3000m"

function makeDirectory {
    local newDirectory=$1

    if !( test -d "${newDirectory}/" ); then
        mkdir "${newDirectory}/"
    fi
}

function unzipGnomad {
    unzipped=$(echo ${gnomad} | sed 's/vcf.bgz/vcf')

    gunzip -c $gnomad \
      > $unzipped
    
    export gnomad=$unzipped
}

function indexFeatureFile {
    gatk --java-options "${javaOpt}" \
      IndexFeatureFile \
        -I ${gnomad}
}

# To take only variants in the interval
function selectVariants {
    makeDirectory $wd

    gatk --java-options "${javaOpt}" \
      SelectVariants \
        -V $gnomad \
        -L $regions \
        -O ${wd}selected.vcf \
        --lenient
}

# gnomad input in mutect2:
# clear ID and QUAL fields and delete all INFO fields other than AF
# note that input must be a plain-text vcf, not a vcf.gz.
# this task re-indexes and compresses the output vcf
function makeAlleleFrequencyOnlyVcf {
    # Save off the header for later:
    grep '^#' ${wd}selected.vcf > ${wd}header &

    # Get all lines in the file except the header:
    # Preserve all fields before INFO, Grab only the AF annotation from the INFO Field
    # replace ID (3rd) and QUAL (6th) columns with '.' (empty):
    grep -v "^#" ${wd}selected.vcf | sed -e 's#\(.*\)\t\(.*\)\t\(.*\)\t\(.*\)\t\(.*\)\t\(.*\)\t\(.*\)\t.*;AF=\([0-9]*\.[e0-9+-]*\).*#\1\t\2\t.\t\4\t\5\t.\t\7\tAF=\8#g' > ${wd}simplified_body &

    # Wait for background processes to finish:
    wait

    # Consolidate files:
    cat ${wd}header ${wd}simplified_body > ${wd}AFonly.vcf

    # # Zip the VCF:
    # bgzip ${wd}simplified.vcf -c > ${wd}AFonly.vcf.gz

    # Index output file:
    gatk --java-options "${javaOpt}" \
      IndexFeatureFile \
        -I ${wd}AFonly.vcf

    # Cleanup:
    rm -f ${wd}header ${wd}body ${wd}simplified_info ${wd}simplified_body ${wd}simplified.vcf ${wd}simplified.vcf.idx
}

# Variants for contamination input in Mutect2
function selectCommonBiallelicSNPs {
    gatk --java-options "${javaOpt}" \
      SelectVariants \
        -V ${wd}AFonly.vcf \
        -select-type SNP -restrict-alleles-to BIALLELIC \
        -select "AF > ${minimumAlleleFrequency}" \
        -O ${wd}variants_for_contamination.vcf \
        --lenient
}

selectVariants
makeAlleleFrequencyOnlyVcf
selectCommonBiallelicSNPs
