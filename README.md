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
Uses *parallelProcessing.sh*.
- pairedFastQsToUnmappedBAM
Converts to unmapped .bam files paired .fastq files in the following format:
```
[sample]_[name]_[library]_[direction: R1 or R2]
```
- validateSam
- samToFastqAndBwaMem
- markDuplicates
- sortAndFixTags
- baseRecalibrator
- applyBqsr

#### somaSNP
Uses *parallelSomaticSNP.sh*.

#### deep
Uses *deepvariant.sh*.

#### anno
Uses *annotation.sh*.

#### 2csv
Uses *goCsv.py*.

#### btil
Uses *bed_to_interval_list.sh*.

#### cvfc
Uses *create_variants_for_contamination.sh*.

#### cnvk
Uses *run_cnvkit.sh*.
