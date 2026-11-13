#!/bin/bash

indir=/kuhpc/work/jychoi/j186c375/MIMULUS/FASTQ/ILLUMINA/20260723_MSv2-PE150_NS-TRAP-15/NS-TRAP-15-Choi
outdir=/kuhpc/work/eeb/n873s703/TRAP_seq_spike_5b/trimmed_files

mkdir -p ${outdir}

for fq in ${indir}/*.fastq.gz
do
    base=$(basename ${fq} .fastq.gz)

    bbduk.sh in=${fq} \
        literal=GACAATCCGTCGAGCAGAGTT \
        k=21 mink=8 ktrim=l \
        out=tmp.fq hdist=3 rcomp=f

    bbduk.sh in=tmp.fq \
        literal=GGTTTCAGGTTTCAGGTTTCAGGTTTGAATTCGG \
        k=31 ktrim=r \
        out=tmp1.fq hdist=3 rcomp=f minlen=1

    bbduk.sh in=tmp1.fq \
        literal=CCGAATTCAAACCTGAAACCTGAAACCTGAAACC \
        k=31 mink=12 ktrim=l \
        out=tmp2.fq hdist=3 rcomp=f

    bbduk.sh in=tmp2.fq \
        literal=AACTCTGCTCGACGGATTGTC \
        k=21 ktrim=r \
        out=${outdir}/${base}_trimmed.fastq.gz \
        hdist=3 rcomp=f minlen=1

    rm -f tmp.fq tmp1.fq tmp2.fq

    echo "${base} finished"

done
