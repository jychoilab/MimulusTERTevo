import csv
from trviz.decomposer import Decomposer
from difflib import SequenceMatcher
import gzip
import csv
import sys
from Bio import SeqIO
import argparse

# example
# fastq_gz_file = "/kuhpc/work/jychoi/j186c375/MIMULUS/FASTQ/ILLUMINA/0251113_MSv2-n-PE150_SK-TRAP-6-Choi/SK-TRAP-6-Choi/LF10_TelPR_S1_L001_R1_001.fastq.gz"
# output_header = test_1
# max_chunk_len
# min_chunk_len


# fastq_gz_file = sys.argv[1]
# # default_motif=sys.argv[2] # number of bases we want to check
# output_header = sys.argv[2] if len(sys.argv) > 2 else "output_parse"


# # max_sequence_length, min and max chunk_len, and allowed_mutations are optional parameters 
# max_sequence_length = int(sys.argv[3]) if len(sys.argv) > 3 else 96
# # set defaults values if not provided
# max_chunk_len = int(sys.argv[4]) if len(sys.argv) > 4 else 10
# min_chunk_len = int(sys.argv[5]) if len(sys.argv) > 5 else 4
# # para allowed_mutations parameter (default 2)
# allowed_mutations = int(sys.argv[6]) if len(sys.argv) > 6 else 2

##### feb 24- adjust that to argsparse w flags

parser = argparse.ArgumentParser(description="Process trap product")
# by default always have
parser.add_argument("fastq_gz_file", help="Input FASTQ.gz file")
# optional - if the name too long leme know
parser.add_argument("--output_name", default="output_parse",help="Output file header (default: output_parse)")
parser.add_argument("--max_seq_length", type=int, default=96,help="Maximum sequence length (default: 96)")
parser.add_argument("--max_chunk_len", type=int, default=10,help="Maximum chunk length (default: 10)")
parser.add_argument("--min_chunk_len", type=int, default=4, help="Minimum chunk length (default: 4)")
parser.add_argument("--max_mutation", type=int, default=2,help="Allowed mutations (default: 2)")

args = parser.parse_args()

# Access variables
fastq_gz_file = args.fastq_gz_file
output_header = args.output_name
max_sequence_length = args.max_seq_length
max_chunk_len = args.max_chunk_len
min_chunk_len = args.min_chunk_len
allowed_mutations = args.max_mutation

include_chunks = True  # Always print all columns to TSV file

# List of motifs to examine, 6bbp and 7bp
default_motifs = ["AACCGA", "AACCCGA"]

def parse_fastq(fastq_file):
    sequences = []

    if fastq_file.endswith(".gz"):
        handle = gzip.open(fastq_file, "rt")
    else:
        handle = open(fastq_file, "r")

    with handle:
        for record in SeqIO.parse(handle, "fastq"):
            sequences.append((record.id, str(record.seq)))
    return sequences


def classify_mutation(a, b,chunk_index, total_chunks):
    """
    Classifies the type of mutation between the motif of interest and the TRAP product.
    inputs:
        a (str): Original sequence (motif).
        b (str): Mutated sequence (TRAP product).
    
    then return list: A list of mutation descriptions (substitution, insertion, or deletion).
    """
    mutations = []
    s = SequenceMatcher(None, a, b)
    
    # Iterate through the opcodes to find mutation types
    # for opcode, a0, a1, b0, b1 in s.get_opcodes():
    #     if opcode == 'replace':
    #         mutations.append(f"Substitution: {a[a0:a1]} <= {b[b0:b1]} at pos {a0}")
    #     elif opcode == 'insert':
    #         mutations.append(f"Deletion: {b[b0:b1]} at pos {a0}")
    #     elif opcode == 'delete':
    #         mutations.append(f"Insertion: {a[a0:a1]} at pos {a0}")

    for opcode, a0, a1, b0, b1 in s.get_opcodes():
        # consider not doing this - so if mutation 
        # if chunk_index == 0 and a0 == 0:  # Mutation at the start of the first chunk
        #     continue
        # if chunk_index == total_chunks - 1 and b1 == len(b):  # Mutation at the end of the last chunk
        #     continue

            
        if opcode == 'replace':
            for i in range(a0, a1):
                if i < len(a) and (b0 + (i - a0)) < len(b):
                    mutations.append(f"{i}_S_{b[b0 + (i - a0)]}{a[i]}")
            #mutations.append(f"{a0}_S_{b[b0:b1]}{a[a0:a1]}")
        elif opcode == 'insert':
            for i in range(b0, b1):
                #if i < len(a) and (b0 + (i - a0)) < len(b):
                mutations.append(f"{i}_D_{b[b0:b1]}")
            #mutations.append(f"{a0}_D_{b[b0:b1]}")
        # elif opcode == 'delete':
        #     for i in range(a0, a1):
        #         #if i < len(a) and (b0 + (i - a0)) < len(b):
        #         mutations.append(f"{i}_D_{a[a0:a1]}")
            #mutations.append(f"{a0}_I_{a[a0:a1]}")
        elif opcode == 'delete':
            for i in range(a0, a1):
                mutations.append(f"{i}_I_{a[i]}")  # Deletion
            # If b is shorter than a, treat the remaining portion of 'a' as deletion
        # if len(b) < len(a):
        #     for i in range(len(b), len(a)):
        #         mutations.append(f"{i}_D_{a[i]}")  # Missing from b (deletion)
    
    return mutations

def analyze_sequences(seq_list, default_motif, max_sequence_length=None):
    """
    Analyzes each sequence in the sequence list, determines motif by C/G content,
    and identifies mutations relative to the chosen motif.
    
    seq_list (list): List of tuples (ID, sequence) to be analyzed.
    default_motif (str): Default motif to use if C and G counts are equal.
    
    return list of rows with sequence analysis results (for TSV output).
    """
    all_rows = []
    trans_table = str.maketrans('ACGT', 'TGCA')
    count_A_motif = default_motif.count("A")
    count_T_motif = default_motif.count("T")

    # determine motif 
    if count_A_motif > count_T_motif: # initial motif is A
        motif_A=default_motif
        motif_T=default_motif.translate(trans_table)[::-1]
    if count_A_motif < count_T_motif: # initial motif is T
        motif_T=default_motif
        motif_A=default_motif.translate(trans_table)[::-1]
    print("motif_T", motif_T)
    print("motif_A",motif_A)

    for seq_id, seq in seq_list:
        # Only analyze sequences that are exactly max_sequence_length bp
        print("max_sequence_length",max_sequence_length)
        if max_sequence_length is not None and len(seq) > max_sequence_length:
            continue  # Skip sequences larger than 96 or any value of max seq len. bp
        seq_to_analyze = seq
        print(f"\nWorking on: {seq_id}")

        # Determine motif based on C and G count for the first 50bp
        count_A = seq_to_analyze.count("A")
        count_T = seq_to_analyze.count("T")

        # determine what motif to use in the sequence
        if count_A > count_T:
            motif = motif_A
        elif count_T > count_A:
            motif = motif_T
        else:
            motif = motif_A

        # if default_motif == 6:  # 6bp
        #     if count_A > count_T:
        #         motif = "AACCGA"
        #     elif count_T > count_A:
        #         motif = "TCGGTT"
        #     else:
        #         motif = "AACCGA"  # If equal CG, use the user-provided motif
        
        # if default_motif == 7:  # 7bp
        #     if count_A > count_T:
        #         motif = "AACCCGA"
        #     elif count_T > count_A:
        #         motif = "TCGGGTT"
        #     else:
        #         motif = "AACCCGA"

        # if default_motif != 6 and default_motif != 7:
        #     print("not support k-mer at this point")
        #     print("running with 6-mer (AACCGA)")
        #     if count_A > count_T:
        #         motif = "AACCGA"
        #     elif count_T > count_A:
        #         motif = "TCGGTT"
        #     else:
        #         motif = "AACCGA"
        
        print(f"Motif chosen for this sequence: {motif} (C={count_A}, G={count_T})")

        # cut chunks using only the first 50bp
        tr_decomposer = Decomposer()
        chunks = tr_decomposer.decompose(seq_to_analyze, motif)

        print("chunks",chunks)

        num_chunks = len(chunks)
        all_chunks = []
        all_mutation_details = []
        filtered_chunk=[]
        problematic=[]

        for i, chunk in enumerate(chunks):
            # pilot: Check if the first or last chunk is too short (less than 3 bp)
            # if (i == 0 and len(chunk) < 3) or (i == num_chunks - 1 and len(chunk) < 3):
            #     continue  # Skip the first and last chunks if they're shorter than 3 bp
            
            # get mutation info
            mutations = classify_mutation(chunk, motif, i, num_chunks)
            # print("chunk", chunk)
            # print("mutation",mutations)
            # print("len mutations", len(mutations))

            #mutation_details = ",".join(mutations) if mutations else "None"
            if not mutations:
                mutation_details = "None"                    # no mutations found
            elif len(mutations) <= allowed_mutations:
                mutation_details = ",".join(mutations)       # keep valid mutations
            else:
                mutation_details = "Filtered"             # too many errors (more than allowed)


            all_chunks.append(chunk)
            all_mutation_details.append(mutation_details)

            # print("mutation_details",mutation_details)
            
            # high quality chunks need to have: longer than min chunk len, shorter than max chunk len
            # number of mutations must be fewer than or equal to allowed_mutations
            if len(chunk) > min_chunk_len and len(chunk) < max_chunk_len and len(mutations) <= allowed_mutations or mutation_details == "None":
                filtered_chunk.append(chunk)
            else:
                problematic.append(chunk)
        
        filtered_mut_details = []
        filtered_mut_count = 0
        for i, chunk in enumerate(all_chunks):
            if chunk in filtered_chunk:
                # print(all_chunks)
                # print(filtered_chunk)
                md = all_mutation_details[i]
                #filtered_mut_details.append(md)
                if md != "Filtered":
                    filtered_mut_details.append(md)
                    if md!="None":
                        filtered_mut_count += len(md.split(","))


        # prepare writing out
        row = {
            "Sequence_ID": seq_id,
            "Sequence": seq,
            "Analyzed_Region": seq_to_analyze,  # Add the region actually analyzed
            "Pattern_check": motif,  
            "Chunks": ",".join(all_chunks),
            "Filtered_chunks": ",".join(filtered_chunk),
            "Num_chunks_filtered": len(filtered_chunk),  # Count of chunks after excluding short ones (after filtered)
            "Num_chunks_raw": len(all_chunks),  # Number of raw chunks before filtering
            # "Mutations_count": sum(len(m.split(',')) for m in all_mutation_details if m != "None"),
            # "Mutation_details": ";".join(all_mutation_details),
            "Mutations_count": filtered_mut_count,          # filtered chunks only
            "Mutation_details": ";".join(filtered_mut_details),  # filtered chunks only
            "problematic_ones": ";".join(problematic) # problematic one
        }

        all_rows.append(row)

    return all_rows


def write_merged_results(output_file, all_results, include_chunks=True):
    """
    Merges results from multiple motifs and writes to a single TSV file.
    Each motif gets its own set of columns with mutation count and filtering info on one line.
    
    inputs:
        output_file (str): Path to the output TSV file.
        all_results (list): List of tuples (motif, rows) for each motif analyzed.
        include_chunks (bool): Whether to include Chunks and Filtered_chunks columns in output.

    """
    merged_rows = []
    
    # Get sequence IDs 
    seq_ids = [row["Sequence_ID"] for row in all_results[0][1]]
    
    for seq_id in seq_ids:
        merged_row = {"Sequence_ID": seq_id}
        
        # Get the sequence 
        first_row = next(r for r in all_results[0][1] if r["Sequence_ID"] == seq_id)
        merged_row["Sequence"] = first_row["Sequence"]
        
        # Add columns for each motif (chunks, filtering info, and mutation count)
        for motif_input, rows in all_results:
            row = next(r for r in rows if r["Sequence_ID"] == seq_id)
            
            # columns for either 6 or 7bp (using suffix format)
            motif_suffix = f"_{len(motif_input)}bp"
            if include_chunks:
                merged_row[f"Chunks{motif_suffix}"] = row["Chunks"]
                merged_row[f"Filtered_chunks{motif_suffix}"] = row["Filtered_chunks"]
                # include problematic chunks in full output
                merged_row[f"problematic_ones{motif_suffix}"] = row.get("problematic_ones", "")
            merged_row[f"Num_chunks_raw{motif_suffix}"] = row["Num_chunks_raw"]
            merged_row[f"Num_chunks_filtered{motif_suffix}"] = row["Num_chunks_filtered"]
            merged_row[f"Mutations_count{motif_suffix}"] = row["Mutations_count"]
            merged_row[f"Mutation_details{motif_suffix}"] = row["Mutation_details"]
        
        # compare 6 and 7 bp based on average mutations per filtered chunk
        chunks_6bp = merged_row.get("Num_chunks_filtered_6bp", 0)
        chunks_7bp = merged_row.get("Num_chunks_filtered_7bp", 0)
        mutations_6bp = merged_row.get("Mutations_count_6bp", float('inf'))
        mutations_7bp = merged_row.get("Mutations_count_7bp", float('inf'))

        # Calculate average mutations per filtered chunk (avoid divide by zero)
        avg_mut_6bp = mutations_6bp / chunks_6bp if chunks_6bp else float('inf')
        avg_mut_7bp = mutations_7bp / chunks_7bp if chunks_7bp else float('inf')

        # If both have no filtered chunks, cannot analyze - say NA
        if chunks_6bp == 0 and chunks_7bp == 0:
            best_motif = "NA"
        # If one has no filtered chunks, the other is better 
        elif chunks_6bp == 0:
            best_motif = "TR2"
        elif chunks_7bp == 0:
            best_motif = "TR1"
        # If both have filtered chunks, compare by average mutation per filtered chunk (lower is better)
        else:
            if avg_mut_6bp < avg_mut_7bp:
                best_motif = "TR1" # 6bp is better
            elif avg_mut_7bp < avg_mut_6bp:
                best_motif = "TR2"
            else:
                best_motif = "same"
        
        merged_row["Class"] = best_motif # merge that in
        
        merged_rows.append(merged_row)
    
    # prep to write out
    fieldnames = ["Sequence_ID", "Sequence"]
    for motif_input, _ in all_results:
        motif_suffix = f"_{len(motif_input)}bp"
        if include_chunks:
            fieldnames.extend([
                f"Chunks{motif_suffix}",
                f"Filtered_chunks{motif_suffix}",
                f"problematic_ones{motif_suffix}"
            ])
        fieldnames.extend([
            f"Num_chunks_raw{motif_suffix}",
            f"Num_chunks_filtered{motif_suffix}",
            f"Mutations_count{motif_suffix}",
            f"Mutation_details{motif_suffix}"
        ])
    fieldnames.append("Class") 
    
    # Write data to TSV file
    with open(output_file, mode="w", newline="") as file:
        writer = csv.DictWriter(file, fieldnames=fieldnames, delimiter='\t')
        writer.writeheader() 
        writer.writerows(merged_rows) 



# Parse the sequences from the fastq.gz file
# but should i add parsing from fastq tho
# seq_list = parse_fastq_gz(fastq_gz_file)
seq_list=parse_fastq(fastq_gz_file)

# seq_list = ["CCGGTTTCGGGTTCGTTTTTGGGTTTCGTTTTTC", "TCGGTTCCGGTTCGGTTTCGGGTTTTCGGTTTCGTGTTT"]

# run analysis for each motif
all_results = []
for motif_input in default_motifs:
    print(f"Analyzing motif: {motif_input}")
    rows = analyze_sequences(seq_list, motif_input, max_sequence_length)
    all_results.append((motif_input, rows))

# Merge results 
write_merged_results(f"{output_header}_results.tsv", all_results, True)
print("done, have a nice day")
