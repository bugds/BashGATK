#!/bin/bash

export picard=/opt/gatk4-data-processing/picard-2.16.0/picard.jar
export regions=/home/bioinfuser/NGS/Reference/intervals/2020_02_02/regions.bed
export outputFolder=/home/bioinfuser/NGS/Reference/intervals/2020_02_02/
export refDict=/home/bioinfuser/NGS/Reference/hg38/hg38.dict

java -jar $picard \
  BedToIntervalList \
    I=$regions \
    O=${outputFolder}regions.interval_list \
    SD=$refDict