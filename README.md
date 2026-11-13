# TRAP sequencing analysis scripts

## Adapter and primer trimming was performed using the script 
01trimming_primeradopter.sh

## The trimmed sequencing reads were then aligned against the Escherichia coli K-12 MG1655 reference genome using BWA-MEM2 as scripted in 
02align_ecoli_single_end.sh

## The mean sequencing coverage of the E. coli spike-in reference genome was calculated for each BAM file using script 
03ecoli_coverage.sh

## Reads that did not align to the E. coli reference genome were subsequently extracted in FASTQ format using 
04extract_non_ecoli_reads.sh

## TRAP product decomposition and counting were performed on the extracted reads using the custom Python script 
05run_main_py.sh

## shell script for running Mver and Mlew data
run_main_py_VER.sh and run_TRAP_LEW.sh

## The normalization of TRAP product counts to E. coli genome coverage, statistic analysis, and visualization of data was performed using 
06_normalized_Plots_TRAP_products.py
