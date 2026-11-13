#!/bin/bash
############################################################
# Single-end alignment of TRAP-seq reads to E. coli genome
# Each FASTQ file (R1 and R2) aligned separately.
# Uses bwa-mem2
############################################################

set -euo pipefail

threads=8

REF=/kuhpc/work/eeb/n873s703/TRAP_seq_spike_5b/reference/ecoli_K12_MG1655.fa

FASTQDIR=/kuhpc/work/eeb/n873s703/TRAP_seq_spike_5b/trimmed_files

OUTDIR=/kuhpc/work/eeb/n873s703/TRAP_seq_spike_5b

BAMDIR=${OUTDIR}/bam_single
STATDIR=${OUTDIR}/stats_single
TMPDIR=${OUTDIR}/tmp_sam_single
LOGDIR=${OUTDIR}/logs_single


mkdir -p "${BAMDIR}" "${STATDIR}" "${TMPDIR}" "${LOGDIR}"


############################################################
# Check bwa-mem2 index
############################################################

if [ ! -s "${REF}.0123" ] || [ ! -s "${REF}.bwt.2bit.64" ]; then
    echo "ERROR: bwa-mem2 index missing"
    echo "Run:"
    echo "bwa-mem2 index ${REF}"
    exit 1
fi


############################################################
# Align FASTQ files
############################################################

FAILED_SAMPLES=()

shopt -s nullglob

FASTQ_FILES=(${FASTQDIR}/*.fastq.gz)


if [ ${#FASTQ_FILES[@]} -eq 0 ]; then
    echo "ERROR: No FASTQ files found in ${FASTQDIR}"
    exit 1
fi


for FASTQ in "${FASTQ_FILES[@]}"
do

    SAMPLE=$(basename "${FASTQ}" .fastq.gz)

    echo "=========================================="
    echo "Processing ${SAMPLE}"
    echo "=========================================="


    SAM="${TMPDIR}/${SAMPLE}.sam"

    BAM="${BAMDIR}/${SAMPLE}.sorted.bam"


    ########################################################
    # bwa-mem2 alignment
    ########################################################

    bwa-mem2 mem \
        -t ${threads} \
        "${REF}" \
        "${FASTQ}" \
        > "${SAM}" \
        2> "${LOGDIR}/${SAMPLE}.bwa.log"


    if [ ! -s "${SAM}" ]; then

        echo "ERROR: Empty SAM for ${SAMPLE}"
        FAILED_SAMPLES+=("${SAMPLE}")
        continue

    fi



    ########################################################
    # Convert SAM -> sorted BAM
    ########################################################

    samtools view \
        -@ ${threads} \
        -b "${SAM}" \
    | samtools sort \
        -@ ${threads} \
        -o "${BAM}"


    if [ ! -s "${BAM}" ]; then

        echo "ERROR: BAM creation failed for ${SAMPLE}"
        FAILED_SAMPLES+=("${SAMPLE}")
        rm -f "${SAM}"
        continue

    fi



    ########################################################
    # BAM indexing
    ########################################################

    samtools index "${BAM}"



    ########################################################
    # Statistics
    ########################################################

    samtools flagstat \
        "${BAM}" \
        > "${STATDIR}/${SAMPLE}.flagstat.txt"


    samtools idxstats \
        "${BAM}" \
        > "${STATDIR}/${SAMPLE}.idxstats.txt"



    rm -f "${SAM}"


    echo "${SAMPLE} finished OK"

done



############################################################
# Summary
############################################################

echo "=========================================="

if [ ${#FAILED_SAMPLES[@]} -eq 0 ]; then

    echo "All samples completed successfully"

else

    echo "Failed samples:"
    printf '%s\n' "${FAILED_SAMPLES[@]}"

fi

echo "=========================================="
