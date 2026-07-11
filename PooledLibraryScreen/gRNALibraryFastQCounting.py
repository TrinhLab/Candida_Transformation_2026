# gRNALibraryFastQCounting
#
# Author: Christopher Cotter
# Manuscript: Cellular Stress Tolerance Governs Genetic Transformability in Recalcitrant Candida Species
#  
# Overview: This script takes an input sequencing FASTQ file and gRNA library csv file and outputs gRNA read counts.
#
# Forward and reverse reads were first paired by aligning ends in pRESTO with the AssemblePairs function in Galaxy
# prior to being input into this script.
#
# Reads are filtered based on the presence of the tRNA sequence or the gRNA scaffold sequence such that partial reads
# are filtered out.
#
# gRNA Library Metadata file should include at minimum a gRNA ID and a gRNA sequence column labeled "Seq_ext". Since 
# we used CRISPR-GRIT gRNAs, our "Seq_ext" sequences contained the 20bp spacer and 5bp immediately upstream of the spacer
# 5' end corresponding to part of the CRISPR-GRIT repair template. This 25bp sequence added extra specificity for counting.
#
# Usage Example:
# python count_gRNAs.py -i input_reads.fastq -m library_metadata.csv -o output_counts.csv

import argparse
import pandas as pd
import sys

def readFastq(filename):
    """Reads FASTQ file and remove the special characters"""
    sequences = []
    qualities = []
    with open(filename) as fh:
        while True:
            fh.readline() # skip name line
            seq = fh.readline().rstrip() # read base sequence
            fh.readline() # skip placeholder line
            qual = fh.readline().rstrip() #base quality line
            if len(seq) == 0:
                break
            sequences.append(seq)
            qualities.append(qual)
    return sequences, qualities

def verify(sequence):
    """This code verifies if a sequence is a DNA or RNA"""
    # set the input sequence
    seq = set(sequence)
     
    # confirm if its elements are equal to the set of valid DNA bases
    # Use a union method to ensure the sequence is verified if it does not contain all the bases
    if seq == {"A", "T", "C", "G"}.union(seq):
        return "DNA"
    elif seq == {"A", "U", "C", "G"}.union(seq):
        return "RNA"
    else:
        return "Invalid sequence"

def rev_comp_st(seq):
    '''This function returns a reverse complement 
    of a DNA or RNA strand'''
    verified = verify(seq)
    if verified == "DNA":
        # complement strand
        seq = seq.replace("A", "t").replace("C", "g").replace("T", "a").replace("G", "c")
        seq = seq.upper()
        # reverse strand
        seq = seq[::-1]
        return seq
 
    elif verified == "RNA":
        # complement strand
        seq = seq.replace("A", "u").replace("C", "g").replace("U", "a").replace("G", "c")
        seq = seq.upper()
        # reverse strand
        seq = seq[::-1]
        return seq
    else:
        return "Invalid sequence"

def main():
    # Set up command line argument parsing
    parser = argparse.ArgumentParser(description="Count gRNAs from a FASTQ file based on a metadata CSV.")
    parser.add_argument("-i", "--input_fastq", required=True, help="Path to the input FASTQ file.")
    parser.add_argument("-m", "--metadata_csv", required=True, help="Path to the gRNA library metadata CSV file.")
    parser.add_argument("-o", "--output_csv", required=True, help="Path for the output CSV file with counts.")
    
    args = parser.parse_args()

    print(f"Reading FASTQ file: {args.input_fastq}")
    try:
        seqs, quals = readFastq(args.input_fastq)
    except FileNotFoundError:
        print(f"Error: Could not find FASTQ file at {args.input_fastq}")
        sys.exit(1)

    print(f"Reading Metadata CSV: {args.metadata_csv}")
    try:
        gRNAs = pd.read_csv(args.metadata_csv)
    except FileNotFoundError:
        print(f"Error: Could not find Metadata CSV file at {args.metadata_csv}")
        sys.exit(1)

    if 'Seq_ext' not in gRNAs.columns:
        print("Error: Metadata CSV must contain a column named 'Seq_ext'")
        sys.exit(1)

    # Search sequences for tRNA and gRNA scaffold
    trna = 'CTTATCTCGTCCA'
    scaffold = 'gttttagagc'
    
    # Search sequences reverse complement in case fastq contains forward and reverse reads
    revtrna = rev_comp_st(trna)
    revscaffold = rev_comp_st(scaffold.upper())

    # Print fastq count statistics
    print("\n--- Initial Sequence Statistics ---")
    count_trna = sum(trna in s for s in seqs)
    print("tRNA:", count_trna)
    count_scaff = sum(scaffold.upper() in s for s in seqs)
    print("Scaffold:", count_scaff)
    count_rev_tRNA = sum(revtrna in s for s in seqs)
    print('Reverse tRNA:', count_rev_tRNA)
    count_rev_scaff = sum(revscaffold in s for s in seqs)
    print('Reverse scaffold:', count_rev_scaff)
    print("Total input sequences:", len(seqs))

    # Filter seqs that don't have the tRNA or gRNA scaffold sequence
    print("\nFiltering sequences...")
    new_seqs = []
    for s in seqs:
        if s.count(revtrna) or s.count(revscaffold):
            new_seqs.append(rev_comp_st(s))
        elif s.count(trna) or s.count(scaffold.upper()):
            new_seqs.append(s)

    print("\n--- Filtered Sequence Statistics ---")
    count_new_trna = sum(trna in s for s in new_seqs)
    count_new_scaffold = sum(scaffold.upper() in s for s in new_seqs)
    print("New tRNA:", count_new_trna)
    print("New scaffold:", count_new_scaffold)
    print("Total filtered sequences:", len(new_seqs))
    print("Invalid sequences:", sum('Invalid sequence' in s for s in new_seqs))

    # Append gRNA counts to new column in dataframe
    print("\nCounting gRNAs (this may take a moment)...")
    gRNA_count = []
    for ind in gRNAs.index:
        gRNA_count.append(sum(gRNAs['Seq_ext'][ind] in s for s in new_seqs))
    
    gRNAs['Counts'] = gRNA_count

    # Write gRNA dataframe to csv
    print(f"\nWriting results to: {args.output_csv}")
    gRNAs.to_csv(args.output_csv, index=False)
    print("Done!")

if __name__ == "__main__":
    main()
