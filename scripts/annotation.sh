set -e
set -o pipefail

function makeDirectory {
    local newDirectory=$1

    if !( test -d "${outputFolder}${newDirectory}/" ); then
        mkdir "${outputFolder}${newDirectory}/"
    fi
}

function decomposeNormalize {
    $bcftools \
        norm -m-both \
        -o ${outputFolder}/${3}/${2}_step1 \
        $1

    $bcftools \
        norm \
        -f $refFasta \
        -o ${outputFolder}/${3}/${2}_step2 \
        ${outputFolder}/${3}/${2}_step1
}

function annotateAnnovar {
    $annovar \
        ${outputFolder}/${2}/${1}_step2 \
        $annovarDb \
        -buildver $buildVer \
        -out ${outputFolder}/${2}/${1}_anno \
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
        -i ${outputFolder}/${2}/${1}_anno.hg38_multianno.vcf \
        -o ${outputFolder}/${2}/${1}_vep.vcf
}

function onlyPASS {
    local file=${outputFolder}/${2}/${1}_vep.vcf

    awk -F '\t' '{if($0 ~ /\#/) print; else if($7 == "PASS") print}' $file > ${file}.pass.vcf
}

if [[ -d "${inputFolder}/mutect2" ]]
then
    for file in ${inputFolder}/mutect2/*/*.filtered.vcf; do
        base=$(basename -- $file)
        folder="anno_soma"
        makeDirectory $folder
        decomposeNormalize $file $base $folder
        annotateAnnovar $base $folder
        annotateVep $base $folder
        onlyPASS $base $folder
    done
fi

if [[ -d "${inputFolder}/deepvariant" ]]
then
    for file in ${inputFolder}/deepvariant/*.vcf.gz; do
        gunzip -k $file
    done
    rm -f ${inputFolder}/deepvariant/*.g.vcf
    for file in ${inputFolder}/deepvariant/*.vcf; do
        base=$(basename -- $file)
        folder="anno_germ"
        makeDirectory $folder
        decomposeNormalize $file $base $folder
        annotateAnnovar $base $folder
        annotateVep $base $folder
        onlyPASS $base $folder
    done
fi
