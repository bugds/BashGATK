#!/bin/bash

set -e
set -o pipefail

java -jar $picard \
  BedToIntervalList \
    I=$regions \
    O=${outputFolder}capture_targets.interval_list \
    SD=$refDict
