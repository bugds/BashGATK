#!/bin/bash

export regions=/home/bioinfuser/NGS/Reference/intervals/2020_02_02/capture_targets.bed
# export outputFolder=/home/bioinfuser/NGS/Reference/intervals/2020_02_02/

java -jar $picard \
  BedToIntervalList \
    I=$regions \
    O=${outputFolder}capture_targets.interval_list \
    SD=$refDict
