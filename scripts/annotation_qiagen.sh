set -e
set -o pipefail

function makeDirectory {
    local newDirectory=$1

    if !( test -d "${outputFolder}${newDirectory}/" ); then
        mkdir "${outputFolder}${newDirectory}/"
    fi
}

function liftoverVcf {
    if !( test -f "${outputFolder}/annotation/${2}_lifted.vcf" ); then
        $gatk LiftoverVcf \
            --INPUT $1 \
            --OUTPUT ${outputFolder}/annotation/${2}_lifted.vcf \
            --REFERENCE_SEQUENCE $refFasta \
            --CHAIN $chain \
            --REJECT ${outputFolder}/annotation/${2}_prerejected.vcf
    fi
}

function checkRejected {
    # Coordinates
    grep "MismatchedRefAllele" ${outputFolder}/annotation/$1_prerejected.vcf > ${outputFolder}/annotation/$1_rejected.vcf
    cat ${outputFolder}/annotation/${1}_rejected.vcf | grep -v ^## | awk -F 'AttemptedLocus=' '{print $2}' | cut -d ";" -f 1 | cut -d "=" -f 2 | tail -n +2 | sed "s/:/\t/" | sed "s/-/\t/" > ${outputFolder}/annotation/${1}_rejected.bed
    cat ${outputFolder}/annotation/${1}_rejected.bed | awk -F '\t' '{ printf("%s\t%s\t%d\n", $1, $2-1, $3) }' > ${outputFolder}/annotation/${1}_rejected_plus.bed
    # Nucleotides
    if [[ $(diff <(cat ${outputFolder}/annotation/${1}_rejected.vcf | grep -v ^## | cut -f 5 | tail -n +2) <($bedtools getfasta -fi $refFasta -bed ${outputFolder}/annotation/${1}_rejected_plus.bed | sed -n 'n;p')) ]]; then
        echo "$1 is bad"
        exit 1
    else
        echo "$1 is ok"
    fi
}
export -f checkRejected

function decomposeNormalize {
    $bcftools \
        norm -m-both \
        -o ${outputFolder}/${folder}/${base}_step1 \
        ${outputFolder}annotation/${1}_lifted.vcf

    $bcftools \
        norm \
        -f $refFasta \
        -o ${outputFolder}/${folder}/${base}_step2 \
        ${outputFolder}/${folder}/${base}_step1
}
export -f decomposeNormalize

function annotateAnnovar {
    $annovar \
        ${outputFolder}annotation/${1}_step2 \
        $annovarDb \
        -thread $parallelJobs\
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
export -f annotateAnnovar

function annotateVep {
    base=$(basename -- $1)
    $vep \
        --offline \
        --everything \
        --fasta $refFasta \
        --vcf \
        --fork $parallelJobs\
        --buffer_size $bufferSize \
        --no_stats \
        --force_overwrite \
        --plugin AlphaMissense,file=$AlphaMissensePath \
        -i ${outputFolder}/${folder}/${base}_anno.hg38_multianno.vcf \
        -o ${outputFolder}/${folder}/${base}_vep.vcf
}
export -f annotateVep

function onlyPASS {
    base=$(basename -- $1)
    local file=${outputFolder}/${folder}/${base}_vep.vcf

    awk -F '\t' '{if($0 ~ /\#/) print; else if($7 == "PASS") print}' $file > ${file}.pass.vcf
}
export -f onlyPASS

# function annotate {
#     filename=$1
#     base=$(basename -- $filename)
#     echo $filename
#     liftoverVcf $file $base
#     checkRejected $base
#     decomposeNormalize $filename
#     annotateAnnovar $filename
#     annotateVep $filename
# }

for file in ${outputFolder}annotation/*.anno.vcf; do
    export folder="annotation"
    base=$(basename -- $file)
    echo $base
    liftoverVcf $file $base
    checkRejected $base
    decomposeNormalize $base
    annotateAnnovar $base
    annotateVep $base
    onlyPASS $base
done


