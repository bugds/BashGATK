# bashgatk.sh

export outputFolder='path_to_output_folder'
export conda='path_to_conda'
export nameSubString='substring_indicating_sample_name_delimited_by_underscores'
export laneSubString='substring_indicating_lane_delimited_by_underscores'
export scriptsDirectory='path_to_scripts_directory_from_the_repository'
export gatk='path_to_gatk'
export fastqc='path_to_fastqc'
export fastp='path_to_fastp'
export samtools='path_to_samtools-1.3.1'
export bedtools='path_to_bedtools'
export picard='path_to_picard-2.16.0'
export bwa='path_to_bwa-0.7.15'
export refFasta='path_to_referencial_genome'
export refDict='path_to_referencial_genome_dictionary'
export dbSnpVcf='path_to_dbsnp_146.hg38.vcf.gz'
export dbSnpVcfIdx='path_to_dbsnp_146.hg38.vcf.gz.tbi'
export millisVcf='path_to_Mills_and_1000G_gold_standard.indels.hg38.vcf.gz'
export millisVcfIdx='path_to_Mills_and_1000G_gold_standard.indels.hg38.vcf.gz.tbi'
export indelsVcf='path_to_known_indels.hg38.vcf.gz'
export indelsVcfIdx='path_to_known_indels.hg38.vcf.gz.tbi'
export regions='path_to_bed-file_with_regions_enriched_by_your_NGS_kit'
export compressionLevel=5
export platform='platform_name'

if [[ $command == 'cvfc' ]]; then
    # create_variants_for_contamination_only_af.sh
    export wd='path_to_directory_with_gnomAD_data'
    export gnomad='path_to_somatic-hg38_af-only-gnomad.hg38.vcf.gz'
    export minimumAlleleFrequency=0.05
    export javaOpt="-Xms3000m"
elif [[ $command == 'proc' ]]; then
    # parallel_processing.sh
    export inputFolder="${outputFolder}/fastq/"
    export bwaVersion="$($bwa 2>&1 | grep -e '^Version' | sed 's/Version: //')"
    export bwaCommandline="$bwa mem -K 100000000 -p -v 3 -t 16 -Y $refFasta"
    export optimal_dup_pixel_distance=100 #In general, a pixel distance of 100 is recommended for data generated using unpatterned flowcells (e.g. HiSeq2500) and a pixel distance of 2500 is recommended for patterned flowcells (e.g. NovaSeq/HiSeq4000).
    export parallelJobs1=5
    export parallelJobs2=1
    export parallelJobs3=1
    export parallelJobs4=1
    export parallelJobs5=1
    export parallelJobs6=1
    export parallelJobs7=1
    export parallelJobs8=1
    export parallelJobs9=1
    export parallelJobs10=1
    export lane=1
elif [[ $command == 'procAmp' ]]; then
    # amplicon_based_processing.sh
    export opticalPixelDistance=100
    export inputFolder="${outputFolder}/fastq/"
    export bwaVersion="$($bwa 2>&1 | grep -e '^Version' | sed 's/Version: //')"
    export bwaCommandline="$bwa mem -K 100000000 -p -v 3 -t 16 -Y $refFasta"
    export parallelJobs=3
    export lane=1
elif [[ $command == 'kumi' ]]; then
    # kapaumi.sh
    export fgbio='path_to_fgbio-2.0.2.jar'
    export inputFolder="${outputFolder}/fastq/"
    export bwaVersion="$($bwa 2>&1 | grep -e '^Version' | sed 's/Version: //')"
    export bwaCommandline1="$bwa mem -t 10 -M $refFasta"
    export bwaCommandline2="$bwa mem -v 3 -t 8 -Y -M $refFasta"
    export optimal_dup_pixel_distance=2500 #In general, a pixel distance of 100 is recommended for data generated using unpatterned flowcells (e.g. HiSeq2500) and a pixel distance of 2500 is recommended for patterned flowcells (e.g. NovaSeq/HiSeq4000).
    export parallelJobs1=5
    export parallelJobs2=1
    export parallelJobs3=1
    export parallelJobs4=1
    export parallelJobs5=1
    export parallelJobs6=1
    export parallelJobs7=1
    export parallelJobs8=1
    export parallelJobs9=1
    export parallelJobs10=1
    export parallelJobs11=1
    export parallelJobs12=1
    export parallelJobs13=1
    export parallelJobs14=1
    export parallelJobs15=1
    export parallelJobs16=1
    export parallelJobs17=1
    export lane=1
elif [[ $command == 'deep' ]]; then
    # deepvariant.sh
    export numCpu=16
    export parallelJobs=5
    export binVersion='latest'
    export refFolder='beginning_of_a_path_to_reference'
    export refFastaPathPart='ending_of_a_path_to_reference'
elif [[ $command == 'mutect2' ]]; then
    # run_mutect2.sh
    export refImg='path_to_hg38.fasta.img'
    export gnomad='path_to_gnomAD_vcf_file'
    export variantsForContamination='path_to_variants_for_contamination.vcf'
    export parallelJobs=5
    export pon='path_to_somatic-hg38_1000g_pon.hg38.vcf'
    export javaOpt="-Xms3000m"
elif [[ $command == 'hapcall' ]]; then
    # run_haplotype_caller.sh
    export parallelJobs=5
    export javaOpt="-Xms3000m"
elif [[ $command == 'anno' ]]; then
    # annotation.sh
    export vep='path_to_vep'
    export bufferSize=500
    export inputFolder="${outputFolder}"
    export annovar='path_to_table_annovar.pl'
    export annovarDb='path_to_annovar_database'
    export bcftools='path_to_bcftools'
    export buildVer=hg38
    export protocol='refGene,ensGene,gnomad211_exome,gnomad312_genome,dbnsfp42a,dbscsnv11,clinvar_latest,cosmic97_coding,cosmic97_noncoding,avsnp150'
    export operation='g,g,f,f,f,f,f,f,f,f'
    export xreffile='path_to_gene_fullxref.txt'
    export parallelJobs=12
    export update_annovar_db='path_to_bugds_update.sh'
elif [[ $command == 'qgen' ]]; then
    # annotation_qiagen.sh
    export chain='path_to_hg19ToHg38.over.chain'
elif [[ $command == '2csv' ]]; then
    # go_csv.py --> run_phylongs.py
    export depth_limit=10
    export af_limit=0.02
    export paf_limit=0.05
    export freq_file='path_to_your_freqs.tsv'
    export clinvar='path_to_hg38_clinvar_latest.txt'
    export constraint='path_to_gnomad.v2.1.1.lof_metrics.by_gene.txt'
    export reference='path_to_gencode_v42_genes.gff3'
    export maskGenes='path_to_mask_genes'
    export annotation='grch38'
    export omim='path_to_genemap2_redacted.tsv'
    export knownScores='path_to_spliceAIscores.tsv'
    export phylongs='path_to_varlist.txt'
elif [[ $command == 'aved' ]]; then
    # average_depth.sh
    export parallelJobs=12
elif [[ $command == 'cnvk' ]]; then
    # run_cnvkit.sh
    export cnvkit='python3 path_to_cnvkit.py'
    export regions='path_to_cnvkit_targets.bed'
    export antitarget='path_to_cnvkit_antitargets.bed'
    export mappable='path_to_access-hg38.bed'
    export refFlat='path_to_hg38.refFlat.txt'
    export AFonly='path_to_gnomAD_AFonly_file'
    export LC_NUMERIC='en_US.UTF-8'
    export javaOpt="-Xms3000m"
elif [[ $command == 'exdh' ]]; then
    # run_cnvkit.sh
    export Rscript='path_to_BashGATK/scripts/run_exome_depth.R'
    export exomeTableFilePath='path_to_hg38_ExomeDepth_exons.bed'
fi
