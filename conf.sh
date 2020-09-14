# bashgatk.sh

export scriptsDirectory=/home/bioinfuser/NGS/Pipelines/scripts
export gatk=/opt/gatk-4.1.5.0/gatk
export samtools=/opt/gatk4-data-processing/samtools-1.3.1/samtools
export picard=/opt/gatk4-data-processing/picard-2.16.0/picard.jar
export bwa=/opt/gatk4-data-processing/bwa-0.7.15/bwa
export refFasta=/home/bioinfuser/NGS/Reference/hg38/hg38.fasta
export refDict=/home/bioinfuser/NGS/Reference/hg38/hg38.dict
export dbSnpVcf=/home/bioinfuser/NGS/Reference/hg38/dbsnp138.vcf
export dbSnpVcfIdx=/home/bioinfuser/NGS/Reference/hg38/dbsnp138.vcf.idx
export millisVcf=/home/bioinfuser/NGS/Reference/hg38/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz
export millisVcfIdx=/home/bioinfuser/NGS/Reference/hg38/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz.tbi
export indelsVcf=/home/bioinfuser/NGS/Reference/hg38/known_indels.hg38.vcf.gz
export indelsVcfIdx=/home/bioinfuser/NGS/Reference/hg38/known_indels.hg38.vcf.gz.tbi
export compressionLevel=5
export platform=ILLUMINA

if [ $command == 'proc' ]; then
    # parallelProcessing.sh
    export fastqc='/opt/FastQC/fastqc'
    export trimmomatic='java -jar /opt/Trimmomatic-0.39/trimmomatic-0.39.jar'
    export trimCommandLine='ILLUMINACLIP:/home/bioinfuser/NGS/Reference/adapters/TruSeq3-PE.fa:2:30:10'
    export inputFolder=${outputFolder}/fastq/
    export regions=/home/bioinfuser/NGS/Reference/intervals/2020_02_02/capture_targets.bed
    export bwaVersion="$($bwa 2>&1 | grep -e '^Version' | sed 's/Version: //')"
    export bwaCommandline="$bwa mem -K 100000000 -p -v 3 -t 16 -Y $refFasta"
    export parallelJobs=3
    export lane=1
elif [ $command == 'procWGS' ]; then
    # processingWGS.sh
    export inputFolder=${outputFolder}/fastq/
    export bwaVersion="$($bwa 2>&1 | grep -e '^Version' | sed 's/Version: //')"
    export bwaCommandline="$bwa mem -K 100000000 -p -v 3 -t 16 -Y $refFasta"
elif [ $command == 'somaSNP' ]; then
    # parallelSomaticSNP.sh
    export regions=/home/bioinfuser/NGS/Reference/intervals/2020_02_02/capture_targets.interval_list
    export refImg=/home/bioinfuser/NGS/Reference/hg38/hg38.fasta.img
    export gnomad=/home/bioinfuser/NGS/Reference/intervals/2020_02_02/additional/AFonly.vcf
    export variantsForContamination=/home/bioinfuser/NGS/Reference/intervals/2020_02_02/additional/variants_for_contamination.vcf
    export parallelJobs=5
    export javaOpt="-Xms3000m"
elif [ $command == 'germSNP' ]; then
    # parallelGermlineSNP.sh
    export regions=/home/bioinfuser/NGS/Reference/intervals/2020_02_02/capture_targets.interval_list
    export parallelJobs=5
    export javaOpt="-Xms3000m"
elif [ $command == 'anno' ]; then
    # annotation.sh
    export vep=/home/bioinfuser/NGS/Software/ensembl-vep/vep
    export fasta=/home/bioinfuser/.vep/homo_sapiens/99_GRCh38/Homo_sapiens.GRCh38.dna.toplevel.fa.gz
    export inputFolder=${outputFolder}/mutect2/
    export annovar=/home/bioinfuser/NGS/Software/annovar/table_annovar.pl
    export annovarDb=/home/bioinfuser/NGS/Software/annovar/humandb/
    export buildVer=hg38
    export protocol='refGene,ensGene,gnomad211_exome,gnomad30_genome,dbnsfp33a,dbscsnv11,clinvar_20190305,cosmic90_coding,cosmic90_noncoding,avsnp150'
    export operation='g,g,f,f,f,f,f,f,f,f'
    export xreffile=/home/bioinfuser/NGS/Software/annovar/example/gene_fullxref.txt
elif [ $command == 'btil' ]; then
    # bed_to_interval_list.sh
    export regions=/home/bioinfuser/NGS/Reference/intervals/2020_02_02/capture_targets.bed
elif [ $command == 'cvfc' ]; then
    # create_variants_for_contamination.sh
    export wd=/home/bioinfuser/NGS/Reference/intervals/2020_02_02/additional/
    export regions=/home/bioinfuser/NGS/Reference/intervals/2020_02_02/capture_targets.bed
    export gnomad=/home/bioinfuser/NGS/Reference/hg38/gnomad/gnomad.exomes.r2.1.1.sites.liftover_grch38.vcf
    export minimumAlleleFrequency=0.05
    export javaOpt="-Xms3000m"
elif [ $command == 'cnvk' ]; then
    # run_cnvkit.sh
    export cnvkit='python3 /home/bioinfuser/NGS/Software/cnvkit/cnvkit.py'
    export regions='/home/bioinfuser/NGS/Reference/intervals/2020_02_02/cnvkit_targets.bed'
    export antitarget='/home/bioinfuser/NGS/Reference/intervals/2020_02_02/cnvkit_antitargets.bed'
    export mappable='/home/bioinfuser/NGS/Reference/hg38/access-hg38.bed'
    export refFlat='/home/bioinfuser/NGS/Reference/hg38/hg38.refFlat.txt'
    export AFonly='/home/bioinfuser/NGS/Reference/intervals/2020_02_02/additional/AFonly.vcf'
    export LC_NUMERIC='en_US.UTF-8'
    export javaOpt="-Xms3000m"
elif [ $command == 'deep' ]; then
    # deepvariant.sh
    export numCpu=16
    export regions=/home/bioinfuser/NGS/Reference/intervals/2020_02_02/capture_targets.bed
    export binVersion='0.10.0'
    export refFolder=/home/bioinfuser/NGS/Reference/
    export refFastaPathPart=hg38/hg38.fasta
    export regionsPathPart=intervals/2020_02_02/cnvkit_targets.bed
fi
