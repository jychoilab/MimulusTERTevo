#!/bin/bash

############################################################
# Run TRAP product decomposition on VER non-E. coli reads
############################################################

set -uo pipefail

PYTHON_SCRIPT=/home/n873s703/work/NS-TRAP-4-Choi/new_bbduk_trim/script/main_2.py

INPUT_DIR=/kuhpc/work/eeb/n873s703/TRAP_seq_spike_5b/non_ecoli_reads

OUTPUT_DIR=/kuhpc/work/eeb/n873s703/TRAP_seq_spike_5b/TRAP_reads_results/VER

mkdir -p ${OUTPUT_DIR}


for FASTQ in ${INPUT_DIR}/ver_*.fastq.gz
do

    SAMPLE=$(basename ${FASTQ} .fastq.gz)

    echo "=========================================="
    echo "Processing ${SAMPLE}"
    echo "=========================================="

    python3 ${PYTHON_SCRIPT} \
        ${FASTQ} \
        --output_name ${OUTPUT_DIR}/${SAMPLE} \
        --max_seq_length 96 \
        --max_chunk_len 10 \
        --min_chunk_len 4 \
        --max_mutation 2

    if [ $? -eq 0 ]; then
        echo "${SAMPLE} finished successfully"
    else
        echo "ERROR processing ${SAMPLE}"
    fi

    echo

done


echo "=========================================="
echo "All VER TRAP reads processed"
echo "Results saved in:"
echo ${OUTPUT_DIR}
echo "=========================================="
