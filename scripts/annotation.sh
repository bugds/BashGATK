set -e
set -o pipefail

function makeDirectory {
    local newDirectory=$1

    if !( test -d "${outputFolder}${newDirectory}/" ); then
        mkdir "${outputFolder}${newDirectory}/"
    fi
}

function annotateAnnovar {
    $annovar \
        $1 \
        $annovarDb \
        -buildver $buildVer \
        -out ${outputFolder}/annotation/${2}_anno \
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
        -i ${outputFolder}/annotation/${1}_anno.hg38_multianno.vcf \
        -o ${outputFolder}/annotation/${1}_vep.vcf
}

function onlyPASS {
    local file=${outputFolder}/annotation/${1}_vep.vcf

    awk -F '\t' '{if($0 ~ /\#/) print; else if($7 == "PASS") print}' $file > ${file}.pass.vcf
}

makeDirectory annotation

for file in ${inputFolder}/*/*.filtered.vcf; do
    base=$(basename -- $file)
    annotateAnnovar $file $base
    annotateVep $base
    onlyPASS $base
done
