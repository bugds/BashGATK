#!/bin/bash

/home/bioinfuser/NGS/Pipelines/processing.sh /home/bioinfuser/NGS/Data/2019_01_14/fastq /home/bioinfuser/NGS/Data/2019_01_14 1
wait
/home/bioinfuser/NGS/Pipelines/processing.sh /home/bioinfuser/NGS/Data/2019_01_16/fastq /home/bioinfuser/NGS/Data/2019_01_16 2
