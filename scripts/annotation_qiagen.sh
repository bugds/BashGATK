set -e
set -o pipefail

function makeDirectory {
    local newDirectory=$1

    if !( test -d "${outputFolder}${newDirectory}/" ); then
        mkdir "${outputFolder}${newDirectory}/"
    fi
}

function liftoverVcf {
    if !( test -f "${outputFolder}/annotation/$2_lifted.vcf" ); then
        $gatk LiftoverVcf \
            --INPUT $1 \
            --OUTPUT ${outputFolder}/annotation/$2_lifted.vcf \
            --REFERENCE_SEQUENCE $refFasta \
            --CHAIN $chain \
            --REJECT ${outputFolder}/annotation/$2_prerejected.vcf
    fi
}

function checkRejected {
    # Coordinates
    grep "MismatchedRefAllele" ${outputFolder}/annotation/$1_prerejected.vcf > ${outputFolder}/annotation/$1_rejected.vcf
    cat ${outputFolder}/annotation/$1_rejected.vcf | grep -v ^## | awk -F 'AttemptedLocus=' '{print $2}' | cut -d ";" -f 1 | cut -d "=" -f 2 | tail -n +2 | sed "s/:/\t/" | sed "s/-/\t/" > ${outputFolder}/annotation/$1_rejected.bed
    cat ${outputFolder}/annotation/$1_rejected.bed | awk -F '\t' '{ printf("%s\t%s\t%d\n", $1, $2-1, $3) }' > ${outputFolder}/annotation/$1_rejected_plus.bed
    # Nucleotides
    if [[ $(diff <(cat ${outputFolder}/annotation/$1_rejected.vcf | grep -v ^## | cut -f 5 | tail -n +2) <($bedtools getfasta -fi $refFasta -bed ${outputFolder}/annotation/$1_rejected_plus.bed | sed -n 'n;p')) ]]; then
        echo "$1 is bad"
        exit 1
    else
        echo "$1 is ok"
    fi
}

function annotateAnnovar {
    $annovar \
        ${outputFolder}annotation/$1_lifted.vcf \
        $annovarDb \
        -buildver $buildVer \
        -out ${outputFolder}annotation/${1}_anno \
        -remove \
        -protocol $protocol \
        -operation $operation \
        -nastring . \
        -polish \
        -xreffile $xreffile \
        -vcfinput
}

function annotateVep {
    $vep \
        --offline \
        --everything \
        --fasta $fasta \
        --vcf \
        -i ${outputFolder}annotation/${1}_anno.hg38_multianno.vcf \
        -o ${outputFolder}annotation/${1}_vep.vcf
}

function onlyPASS {
    local file=${outputFolder}annotation/${1}_vep.vcf

    awk -F '\t' '{if($0 ~ /\#/) print; else if($7 == "PASS") print}' $file > ${file}.pass.vcf
}

for file in ${outputFolder}annotation/*.anno.vcf; do
    base=$(basename -- $file)
    echo $base
    liftoverVcf $file $base
    checkRejected $base
    annotateAnnovar $base
    annotateVep $base
    onlyPASS $base
done


