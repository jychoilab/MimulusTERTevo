#!/bin/bash

BAMDIR=/kuhpc/work/eeb/n873s703/TRAP_seq_spike_5b/bam_single
OUTDIR=/kuhpc/work/eeb/n873s703/TRAP_seq_spike_5b/ecoli_coverage

mkdir -p ${OUTDIR}

echo -e "Sample\tMean_Ecoli_Coverage" > ${OUTDIR}/ecoli_coverage.tsv


for BAM in ${BAMDIR}/*.sorted.bam
do

    SAMPLE=$(basename ${BAM} .sorted.bam)

    echo "Processing ${SAMPLE}"

    COVERAGE=$(samtools depth -aa ${BAM} | \
    awk '{sum+=$3; count+=1} END {print sum/count}')

    echo -e "${SAMPLE}\t${COVERAGE}" >> ${OUTDIR}/ecoli_coverage.tsv

done
