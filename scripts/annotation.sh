export vep=/home/bioinfuser/NGS/Software/ensembl-vep/vep
export fasta=/home/bioinfuser/.vep/homo_sapiens/99_GRCh38/Homo_sapiens.GRCh38.dna.toplevel.fa.gz
export inputFolder=${outputFolder}/mutect2/
export annovar=/home/bioinfuser/NGS/Software/annovar/table_annovar.pl
export annovarDb=/home/bioinfuser/NGS/Software/annovar/humandb/
export buildVer=hg38
export protocol='refGene,ensGene,gnomad211_exome,gnomad30_genome,dbnsfp33a,dbscsnv11,clinvar_20190305,cosmic90_coding,cosmic90_noncoding,avsnp150'
export operation='g,g,f,f,f,f,f,f,f,f'
export xreffile=/home/bioinfuser/NGS/Software/annovar/example/gene_fullxref.txt

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

makeDirectory annotation

for file in ${inputFolder}/*/*.filtered.vcf; do
    base=$(basename -- $file)
    annotateAnnovar $file $base
    annotateVep $base
done
