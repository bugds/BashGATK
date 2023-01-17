set -e
set -o pipefail

function parallelRun {  
    local files=$2

    export -f $1
    parallel -j $parallelJobs $1 ::: $files
}

function makeDirectory {
    local newDirectory=$1

    if !( test -d "${outputFolder}${newDirectory}/" ); then
        mkdir "${outputFolder}${newDirectory}/"
    fi
}
export -f makeDirectory

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
export -f decomposeNormalize

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
export -f annotateAnnovar

function annotateVep {
    $vep \
        --offline \
        --everything \
        --fasta $fasta \
        --vcf \
        -i ${outputFolder}/${2}/${1}_anno.hg38_multianno.vcf \
        -o ${outputFolder}/${2}/${1}_vep.vcf
}
export -f annotateVep

function onlyPASS {
    local file=${outputFolder}/${2}/${1}_vep.vcf

    awk -F '\t' '{if($0 ~ /\#/) print; else if($7 == "PASS") print}' $file > ${file}.pass.vcf
}
export -f onlyPASS

function annotate {
    base=$(basename -- $1)
    makeDirectory $folder
    decomposeNormalize $1 $base $folder
    annotateAnnovar $base $folder
    annotateVep $base $folder
    onlyPASS $base $folder
}

bash ${update_annovar_db}

if [[ -d "${inputFolder}/mutect2" ]]
then
    folder="anno_soma"
    parallelRun annotate ${inputFolder}/mutect2/*/*.filtered.vcf
fi

if [[ -d "${inputFolder}/deepvariant" ]]
then
    folder="anno_germ"
    gunzip -k -f ${inputFolder}/deepvariant/*.vcf.gz
    rm -f ${inputFolder}/deepvariant/*.g.vcf
    parallelRun annotate ${inputFolder}/mutect2/*/*.filtered.vcf
fi
