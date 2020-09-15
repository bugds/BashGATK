# BashGATK
Bash implementation of GATK-based analysis, designed specifically for paired targeted NGS.

Basic somatic pipeline (represented by *main.sh*):
- bashgatk.sh proc [output folder]
- bashgatk.sh somaSNP [output folder]
- bashgatk.sh anno [output folder]
- bashgatk.sh 2csv [output folder]

bashgatk.sh consists of instructions for each BashGATK tool and general parameters, which are *export*ed to the environment.

Usage:
```
bashgatk.sh [tool name] [output folder]
```

Dependables:
- GATK;
- Samtools;
- BWA;
- Picard.

## Tools

#### proc
Performs preprocessing steps for further analysis.
- pairedFastQsToUnmappedBAM
Converts to unmapped .bam files paired .fastq files in the following format:
```
[sample]_[name]_[library]_[direction: R1 or R2]
```
- validateSam
Checks if the resulting .bam file is good for GATK analysis
- samToFastqAndBwaMem
- markDuplicates
- sortAndFixTags
- baseRecalibrator
- applyBqsr

#### somaSNP
Uses *parallelSomaticSNP.sh*.
- mutect2
- learnReadOrientationModel
- calculateContamination
- filterMutectCalls

#### deep
Uses *deepvariant.sh*.
Runs deepvariant in a Docker container

#### anno
Uses *annotation.sh*.
- ANNOVAR
- VEP
- Filtering to get only *PASS*ed variants

#### 2csv
Uses *goCsv.py*.
Converts vcf to human-readable csv format.

#### btil
Uses *bed_to_interval_list.sh*.
Converts .bed file (capture targets) to interval list.

#### cvfc
Uses *create_variants_for_contamination.sh*.
Given a gnomAD .vcf file, produces files for mutect2:
- sites-only .vcf file where the only *INFO* field is allele frequency (AF),
- the same file restricted to a given minimum allele frequency.

#### cnvk
Uses *run_cnvkit.sh*.
Runs CNVKit to derive CNV information from targeted NGS run.
