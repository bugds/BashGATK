export vep='/home/NGS/Software/ensembl-vep/vep'
export fasta='/home/.vep/homo_sapiens/99_GRCh38/Homo_sapiens.GRCh38.dna.toplevel.fa.gz'
export outputFolder="$(realpath $1)/"

# function annotateAnnovar {

# Gene-based
# refseq
# ensembl

# Filter-based

# frequencies:
# 1000g2015aug
# kaviar_20150923
# hrcr1
# gnomad_genome
# exac03
# esp6500siv2_all
# gnomad_exome

# predictions:
# dbnsfp33a (exome)
# dbscsnv11 (splice)
# spidex (splice)
# genome predictions are too big

# disease-specific:
# clinvar_20190305
# cosmic70
# cosmic90
# nci60

# variant identifiers:
# avsnp142
# }

function annotateVep {

    for file in $files; do
        $vep \
            --cache \
            --everything \
            --fasta $fasta \
            -i ~/NGS/Data/2019_01_15/MED06.filtered.vcf \
            -o ~/Downloads/output.vcf
    done
}