#!/bin/bash

set -e
set -o pipefail

export command=$1
export config=$2

source $config

if [ $command == 'btil' ]; then
    bash ${scriptsDirectory}/bed_to_interval_list.sh
elif [ $command == 'cvfc' ]; then
    bash ${scriptsDirectory}/create_variants_for_contamination_only_af.sh
elif [ $command == 'proc' ]; then
    source $conda
    conda activate fastqc
    bash ${scriptsDirectory}/parallel_processing.sh
elif [ $command == 'procAmp' ]; then
    source $conda
    conda activate fastqc
    bash ${scriptsDirectory}/amplicon_based_processing.sh
elif [ $command == 'kumi' ]; then
    source $conda
    conda activate fastqc
    bash ${scriptsDirectory}/kapaumi.sh
elif [ $command == 'deep' ]; then
    bash ${scriptsDirectory}/deepvariant.sh
elif [ $command == 'mutect2' ]; then
    bash ${scriptsDirectory}/run_mutect2.sh
elif [ $command == 'hapcall' ]; then
    bash ${scriptsDirectory}/run_haplotype_caller.sh
elif [ $command == 'anno' ]; then
    bash ${scriptsDirectory}/annotation.sh
elif [ $command == 'qgen' ]; then
    bash ${scriptsDirectory}/annotation_qiagen.sh
elif [ $command == '2csv' ]; then
    source $conda
    conda activate stat
    python3 ${scriptsDirectory}/go_csv.py
    python3 ${scriptsDirectory}/2rus.py
    python3 ${scriptsDirectory}/create_gene_mask.py
    conda activate spliceai
    python3 ${scriptsDirectory}/run_spliceai.py
    conda activate stat
    python3 ${scriptsDirectory}/add_omim.py
    python3 ${scriptsDirectory}/run_phylongs.py
elif [ $command == '2csvq' ]; then
    source $conda
    conda activate stat
    python3 ${scriptsDirectory}/go_csv_qgen.py
elif [ $command == 'aved' ]; then
    bash ${scriptsDirectory}/average_depth.sh
elif [ $command == 'cnvk' ]; then
    bash ${scriptsDirectory}/run_cnvkit.sh
fi