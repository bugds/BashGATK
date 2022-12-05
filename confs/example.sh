# bashgatk.sh

export scriptsDirectory=/home/bioinfuser/applications/BashGATK/scripts
export gatk=/home/bioinfuser/applications/gatk-4.2.5.0/gatk
export samtools=/home/bioinfuser/applications/samtools-1.3.1/samtools
export picard=/home/bioinfuser/applications/picard-2.16.0/picard.jar
export bwa=/home/bioinfuser/applications/bwa-0.7.15/bwa
export refFasta=/home/bioinfuser/data/hg38/hg38.fasta
export refDict=/home/bioinfuser/data/hg38/hg38.dict
export dbSnpVcf=/home/bioinfuser/data/hg38/dbsnp138.vcf
export dbSnpVcfIdx=/home/bioinfuser/data/hg38/dbsnp138.vcf.idx
export millisVcf=/home/bioinfuser/data/hg38/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz
export millisVcfIdx=/home/bioinfuser/data/hg38/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz.tbi
export indelsVcf=/home/bioinfuser/data/hg38/known_indels.hg38.vcf.gz
export indelsVcfIdx=/home/bioinfuser/data/hg38/known_indels.hg38.vcf.gz.tbi
export compressionLevel=5
export platform=ILLUMINA

if [[ $command == 'proc' ]]; then
    # parallelProcessing.sh
    export opticalPixelDistance=100 #In general, a pixel distance of 100 is recommended for data generated using unpatterned flowcells (e.g. HiSeq2500) and a pixel distance of 2500 is recommended for patterned flowcells (e.g. NovaSeq/HiSeq4000).
    export nameSubString='N'
    export laneSubString='L'
    export fastqc='/home/bioinfuser/applications/FastQC/fastqc'
    export trimmomatic='java -jar /home/bioinfuser/applications/Trimmomatic-0.39/trimmomatic-0.39.jar'
    export trimCommandLine='CROP:20 HEADCROP:20 ILLUMINACLIP:/home/bioinfuser/applications/Trimmomatic-0.39/adapters/TruSeq3-PE-3.fa:2:30:10'
    export inputFolder=${outputFolder}/fastq/
    export regions=/home/bioinfuser/data/brca_seq/intervals/BRCACNVbr283_ROI.bed
    export bwaVersion="$($bwa 2>&1 | grep -e '^Version' | sed 's/Version: //')"
    export bwaCommandline="$bwa mem -K 100000000 -p -v 3 -t 16 -Y $refFasta"
    export parallelJobs=3
    export lane=1
elif [[ $command == 'procAmp' ]]; then
    # parallelProcessingAmpliconBased.sh
    export opticalPixelDistance=100
    export nameSubString='N'
    export laneSubString='L'
    export fastqc='/home/bioinfuser/applications/FastQC/fastqc'
    export trimmomatic='java -jar /home/bioinfuser/applications/Trimmomatic-0.39/trimmomatic-0.39.jar'
    export trimCommandLine='ILLUMINACLIP:/home/bioinfuser/applications/Trimmomatic-0.39/adapters/TruSeq3-PE-3.fa:2:30:10 LEADING:30 TRAILING:30 SLIDINGWINDOW:20:35 MINLEN:120'
    export inputFolder=${outputFolder}/fastq/
    export regions=/home/bioinfuser/data/brca_seq/intervals/BRCACNVbr283_ROI.bed
    export bwaVersion="$($bwa 2>&1 | grep -e '^Version' | sed 's/Version: //')"
    export bwaCommandline="$bwa mem -K 100000000 -p -v 3 -t 16 -Y $refFasta"
    export parallelJobs=3
    export lane=1
elif [[ $command == 'procWGS' ]]; then
    # processingWGS.sh
    export inputFolder=${outputFolder}/fastq/
    export bwaVersion="$($bwa 2>&1 | grep -e '^Version' | sed 's/Version: //')"
    export bwaCommandline="$bwa mem -K 100000000 -p -v 3 -t 16 -Y $refFasta"
elif [[ $command == 'somaSNP' ]]; then
    # parallelSomaticSNP.sh
    export regions=/home/bioinfuser/data/bykova_more/capture_targets.interval_list
    export refImg=/home/bioinfuser/data/hg38/hg38.fasta.img
    export gnomad=/home/bioinfuser/data/bykova_more/AFonly.vcf
    export variantsForContamination=/home/bioinfuser/data/bykova_more/variants_for_contamination.vcf
    export parallelJobs=5
    export pon=/home/bioinfuser/data/hg38/somatic-hg38_1000g_pon.hg38.vcf
    export javaOpt="-Xms3000m"
elif [[ $command == 'germSNP' ]]; then
    # parallelGermlineSNP.sh
    export regions=/home/bioinfuser/data/brca_seq/intervals/capture_targets.interval_list
    export parallelJobs=5
    export javaOpt="-Xms3000m"
elif [[ $command == 'anno' ]]; then
    # annotation.sh
    export vep=/home/bioinfuser/applications/ensembl-vep/vep
    export fasta=/home/bioinfuser/.vep/homo_sapiens/101_GRCh38/Homo_sapiens.GRCh38.dna.primary_assembly.fa
    export inputFolder=${outputFolder}
    export annovar=/home/bioinfuser/applications/annovar/table_annovar.pl
    export annovarDb=/home/bioinfuser/applications/annovar/humandb/
    export bcftools=/home/bioinfuser/applications/bcftools-1.15/bin/bcftools
    export buildVer=hg38
    export protocol='refGene,ensGene,gnomad211_exome,gnomad30_genome,dbnsfp42a,dbscsnv11,clinvar_latest,cosmic95_coding,cosmic95_noncoding,avsnp150'
    export operation='g,g,f,f,f,f,f,f,f,f'
    export xreffile=/home/bioinfuser/applications/annovar/example/gene_fullxref.txt
elif [[ $command == 'btil' ]]; then
    # bed_to_interval_list.sh
    export regions=/home/bioinfuser/data/brca_seq/intervals/BRCACNVbr283_ROI.bed
elif [[ $command == 'cvfc' ]]; then
    # create_variants_for_contamination.sh
    export wd=/home/bioinfuser/data/brca_seq/intervals
    export regions=/home/bioinfuser/data/brca_seq/intervals/BRCACNVbr283_ROI.bed
    export gnomad=/home/bioinfuser/data/hg38/gnomad/gnomad.exomes.r2.1.1.sites.liftover_grch38.vcf
    export minimumAlleleFrequency=0.05
    export javaOpt="-Xms3000m"
elif [[ $command == 'cnvk' ]]; then
    # run_cnvkit.sh
    export cnvkit='python3 /home/bioinfuser/applications/cnvkit/cnvkit.py'
    export regions='/home/bioinfuser/data/MED/intervals/cnvkit_targets.bed'
    export antitarget='/home/bioinfuser/data/MED/intervals/cnvkit_antitargets.bed'
    export mappable='/home/bioinfuser/data/MED/intervals/access-hg38.bed'
    export refFlat='/home/bioinfuser/data/MED/intervals/hg38.refFlat.txt'
    export AFonly='/home/bioinfuser/data/MED/intervals/additional/AFonly.vcf'
    export LC_NUMERIC='en_US.UTF-8'
    export javaOpt="-Xms3000m"
elif [[ $command == 'deep' ]]; then
    # deepvariant.sh
    export numCpu=16
    export regions=/home/bioinfuser/data/brca_seq/intervals/BRCACNVbr283_ROI.bed
    export binVersion='latest'
    export refFolder=/home/bioinfuser/
    export refFastaPathPart=data/hg38/hg38.fasta
    export regionsPathPart=data/brca_seq/intervals/BRCACNVbr283_ROI.bed
elif [[ $command == 'cint' ]]; then
    # createIntervals.sh
    export bedtools=/home/bioinfuser/applications/bedtools.static.binary
elif [[ $command == 'aved' ]]; then
    # averageDepth.sh
    export regions=/home/bioinfuser/data/brca_seq/intervals/BRCACNVbr283_ROI.bed
elif [[ $command == 'qgen' ]]; then
    # annotationQiagen.sh
    export chain=/home/bioinfuser/data/hg38/hg19ToHg38.over.chain
    export bedtools=/home/bioinfuser/applications/bedtools.static.binary
fi
