#!/bin/bash

set -e
set -o pipefail

export command=$1
export config=$2

source $config

if [ $command == 'proc' ]; then
    bash ${scriptsDirectory}/parallelProcessing.sh
elif [ $command == 'procAmp' ]; then
    bash ${scriptsDirectory}/parallelProcessingAmpliconBased.sh
elif [ $command == 'procWGS' ]; then
    bash ${scriptsDirectory}/processingWGS.sh
elif [ $command == 'somaSNP' ]; then
    bash ${scriptsDirectory}/parallelSomaticSNP.sh
elif [ $command == 'germSNP' ]; then
    bash ${scriptsDirectory}/parallelGermlineSNP.sh
elif [ $command == 'anno' ]; then
    bash ${scriptsDirectory}/annotation.sh
elif [ $command == 'qgen' ]; then
    bash ${scriptsDirectory}/annotation_qiagen.sh
elif [ $command == 'btil' ]; then
    bash ${scriptsDirectory}/bed_to_interval_list.sh
elif [ $command == 'cvfc' ]; then
    bash ${scriptsDirectory}/create_variants_for_contamination.sh
elif [ $command == '2csv' ]; then
    source ~/applications/miniconda3/etc/profile.d/conda.sh
    conda activate stat
    python3 ${scriptsDirectory}/goCsv.py
elif [ $command == '2csvq' ]; then
    source ~/applications/miniconda3/etc/profile.d/conda.sh
    conda activate stat
    python3 ${scriptsDirectory}/goCsvQiagen.py
elif [ $command == 'cnvk' ]; then
    bash ${scriptsDirectory}/run_cnvkit.sh
elif [ $command == 'deep' ]; then
    bash ${scriptsDirectory}/deepvariant.sh
elif [ $command == 'deep19' ]; then
    bash ${scriptsDirectory}/deepvariant.sh
elif [ $command == 'cnvk19' ]; then
    bash ${scriptsDirectory}/run_cnvkit.sh
elif [ $command == '2rus' ]; then
    source ~/applications/miniconda3/etc/profile.d/conda.sh
    conda activate stat
    python3 ${scriptsDirectory}/2rus.py
elif [ $command == 'cint' ]; then
    bash ${scriptsDirectory}/createIntervals.sh
elif [ $command == 'aved' ]; then
    bash ${scriptsDirectory}/averageDepth.sh
elif [ $command == 'afcvfc' ]; then
    bash ${scriptsDirectory}/create_variants_for_contamination_only_af.sh $outputFolder
elif [ $command == 'kumi' ]; then
    bash ${scriptsDirectory}/kapaumi.sh
# elif [ $command == 'acmg' ]; then
#     source ~/applications/miniconda3/etc/profile.d/conda.sh
#     conda activate stat
#     python3 ${scriptsDirectory}/acmg_classifier.py
elif [ $command == 'mask' ]; then
    source ~/applications/miniconda3/etc/profile.d/conda.sh
    conda activate stat
    python3 ${scriptsDirectory}/create_gene_mask.py
elif [ $command == 'spai' ]; then
    source ~/applications/miniconda3/etc/profile.d/conda.sh
    conda activate spliceai
    python3 ${scriptsDirectory}/run_spliceai.py
fi
