set -e
set -o pipefail

function makeDirectory {
    local newDirectory=$1

    if !( test -d "${outputFolder}${newDirectory}/" ); then
        mkdir "${outputFolder}${newDirectory}/"
    fi
}
export -f makeDirectory

function decomposeNormalize {
    base=$(basename -- $1)
    $bcftools \
        norm -m-both \
        -o ${outputFolder}/${folder}/${base}_step1 \
        $1

    $bcftools \
        norm \
        -f $refFasta \
        -o ${outputFolder}/${folder}/${base}_step2 \
        ${outputFolder}/${folder}/${base}_step1
}
export -f decomposeNormalize

function annotateAnnovar {
    base=$(basename -- $1)
    $annovar \
        ${outputFolder}/${folder}/${base}_step2 \
        $annovarDb \
        -thread $parallelJobs\
        -buildver $buildVer \
        -out ${outputFolder}/${folder}/${base}_anno \
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

function annotate {
    filename=$1
    echo $filename
    decomposeNormalize $filename
    annotateAnnovar $filename
    annotateVep $filename
}

bash ${update_annovar_db}

echo "Running with ${parallelJobs} parallel threads"

if [[ -d "${outputFolder}/mutect2" ]]
then
    export folder="anno_soma"
    makeDirectory $folder
    files=${outputFolder}/mutect2/*/*.filtered.vcf
    for filename in $files; do
        annotate $filename
        onlyPASS $filename
    done
    rm ${outputFolder}/${folder}/*step1
    rm ${outputFolder}/${folder}/*step2
    rm ${outputFolder}/${folder}/*anno.vcf
    rm ${outputFolder}/${folder}/*anno.txt
    rm ${outputFolder}/${folder}/*anno.avinput
fi

if [[ -d "${outputFolder}/deepvariant" ]]
then
    export folder="anno_germ"
    makeDirectory $folder
    gunzip -k -f ${outputFolder}/deepvariant/*.vcf.gz
    rm -f ${outputFolder}/deepvariant/*.g.vcf
    files=${outputFolder}/deepvariant/*.vcf
    for filename in $files; do
        annotate $filename
        onlyPASS $filename
    done
    rm ${outputFolder}/${folder}/*step1
    rm ${outputFolder}/${folder}/*step2
    rm ${outputFolder}/${folder}/*anno.vcf
    rm ${outputFolder}/${folder}/*anno.txt
    rm ${outputFolder}/${folder}/*anno.avinput
fi

if [[ -d "${outputFolder}/haplotype_caller" ]]
then
    export folder="anno_hapc"
    makeDirectory $folder
    files=${outputFolder}/haplotype_caller/*.vcf
    for filename in $files; do
        annotate $filename
    done
    rm ${outputFolder}/${folder}/*step1
    rm ${outputFolder}/${folder}/*step2
    rm ${outputFolder}/${folder}/*anno.vcf
    rm ${outputFolder}/${folder}/*anno.txt
    rm ${outputFolder}/${folder}/*anno.avinput
fi