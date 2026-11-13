#!/bin/bash

BAMDIR=/kuhpc/work/eeb/n873s703/TRAP_seq_spike_5b/bam_single
OUTDIR=/kuhpc/work/eeb/n873s703/TRAP_seq_spike_5b/non_ecoli_reads

mkdir -p ${OUTDIR}

for BAM in ${BAMDIR}/*.sorted.bam
do

    SAMPLE=$(basename ${BAM} .sorted.bam)

    echo "======================================"
    echo "Extracting unmapped reads: ${SAMPLE}"
    echo "======================================"

    samtools fastq \
        -f 4 \
        ${BAM} \
        | gzip > ${OUTDIR}/${SAMPLE}_unmapped.fastq.gz

    echo "${SAMPLE} finished"

done

echo "All unmapped reads extracted."
