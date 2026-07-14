# CRISPR-GRIT Guide Generator Script
#
# Author: Christopher Cotter
# Manuscript: Cellular Stress Tolerance Governs Genetic Transformability in Recalcitrant Candida Species
#
# Overview: This script reads an input FASTA file and a CSV of gRNAs, aligns the guides, 
# selects homology arms, introduces in-frame stop mutations combined with a frameshift, 
# and outputs the CRISPR-GRIT fragment with BsmBI sites and a 20bp seed.
#
# Usage Example:
# python CRISPRGRIT_GuideGenerator_v1.py -i guides.csv -f gene.fasta -o output.csv --up_overhang GTCCA --down_overhang GTTTT --ha_length 50 --overlap_length 20

import argparse
import sys
import pandas as pd
import numpy as np
from math import ceil
from Bio import SeqIO, pairwise2
from Bio.Seq import MutableSeq, Seq
from Bio.SeqUtils import MeltingTemp as mt

def process_sequence(record, sequence):
    """Aligns the sequence to the FASTA record and returns alignment details."""
    sequence = Seq(sequence)

    for is_direct in [True, False]:
        seq_to_align = sequence if is_direct else sequence.reverse_complement()
        # Strict alignment for CRISPR-GRIT (requires perfect 20bp match)
        aln = pairwise2.align.localms(record.seq, seq_to_align, 5, -4, -2, -0.5)

        if aln:
            best_alignment = aln[0]
            score = best_alignment[2]
            if score == 100:  # 5 * 20bp = 100
                return best_alignment, is_direct
    return None, False

def main():
    parser = argparse.ArgumentParser(description="Generate CRISPR-GRIT guides with homology arms and in-frame stop/frameshift mutations.")
    
    # Required Arguments
    parser.add_argument("-i", "--input_csv", required=True, help="Path to input gRNA CSV file (must contain a 'Sequence' column).")
    parser.add_argument("-f", "--fasta", required=True, help="Path to input FASTA file.")
    parser.add_argument("-o", "--output_csv", required=True, help="Path for the output CSV file.")
    parser.add_argument("--up_overhang", required=True, help="5bp Upstream Golden Gate Overhang (e.g., GTCCA).")
    parser.add_argument("--down_overhang", required=True, help="5bp Downstream Golden Gate Overhang (e.g., GTTTT).")
    parser.add_argument("--ha_length", type=int, required=True, help="Desired homology arm length (e.g., 50).")
    parser.add_argument("--overlap_length", type=int, required=True, help="Desired oligo overlap length (e.g., 20).")
    
    # Optional Arguments
    parser.add_argument("--orf_start", type=int, default=1000, help="Coordinate where the ORF starts in the FASTA (default: 1000).")
    parser.add_argument("--insert_seq", type=str, default="taataa", help="Sequence to insert (default: 'taataa').")

    args = parser.parse_args()

    # Load Data
    print(f"Reading gRNA CSV: {args.input_csv}")
    try:
        df = pd.read_csv(args.input_csv, delimiter=',')
    except FileNotFoundError:
        print(f"Error: Could not find CSV file at {args.input_csv}")
        sys.exit(1)

    if 'Sequence' not in df.columns:
        print("Error: Input CSV must contain a column named 'Sequence'")
        sys.exit(1)

    print(f"Reading FASTA file: {args.fasta}")
    try:
        list(SeqIO.parse(args.fasta, 'fasta'))
    except FileNotFoundError:
        print(f"Error: Could not find FASTA file at {args.fasta}")
        sys.exit(1)

    # Initialize Lists
    GRIT_List = []
    FW_Primer_List = []
    RV_Primer_List = []
    Tm_List = []
    BindingRegion = []
    
    stoptag = args.insert_seq
    BsmBI_Up = 'AAACGTCTCT' + args.up_overhang
    BsmBI_Down = args.down_overhang + 'GGAGACGAAA'
    HA_length = args.ha_length
    primer_len = args.overlap_length
    my_orf_coords = args.orf_start

    print("Processing sequences...")
    for sequence in df['Sequence']:
        found_in_any_record = False
        
        for record in SeqIO.parse(args.fasta, 'fasta'):
            # If the user didn't specify an ORF start, try to find the first ATG
            if args.orf_start == 1000: 
                atg_idx = record.seq.find('ATG')
                if atg_idx != -1:
                    my_orf_coords = atg_idx

            aln_result, is_direct = process_sequence(record, sequence)
            
            if aln_result:
                align_start = aln_result[3]
                align_end = aln_result[4]
                seqstart = align_start
                subseq = record.seq[align_start:align_end]
                
                found = False
                
                if is_direct:
                    CutSiteCoords = seqstart + len(subseq) + 6
                    AA_low = (seqstart + len(subseq) - 4 - my_orf_coords) // 3
                    AA_high = ceil((CutSiteCoords - my_orf_coords) / 3)
                    
                    for i in range(1, 4):
                        Mut = MutableSeq(record.seq)
                        Mut[seqstart + len(subseq) - (i + 1):CutSiteCoords - (i - 1)] = stoptag
                        Mut_AA = Mut[my_orf_coords:].translate(table=12)
                        
                        if Mut_AA[AA_low:AA_high].find('**') != -1:
                            GRIT = Seq(BsmBI_Up + Mut[seqstart + len(subseq) - (i + 1) - HA_length:CutSiteCoords - (i + 1) + HA_length] + sequence + BsmBI_Down)
                            GRIT_List.append(GRIT)
                            
                            mid = len(GRIT) // 2
                            FW_Primer_List.append(Seq(GRIT[:mid + primer_len // 2]))
                            RV_Primer_List.append(Seq(GRIT[mid - primer_len // 2:]).reverse_complement())
                            Tm_List.append(mt.Tm_Wallace(GRIT[mid - primer_len // 2:mid + primer_len // 2]))
                            BindingRegion.append(GRIT[mid - primer_len // 2:mid + primer_len // 2])
                            
                            found = True
                            found_in_any_record = True
                            break
                            
                else:
                    CutSiteCoords = seqstart - 6
                    AA_low = (CutSiteCoords - my_orf_coords) // 3
                    AA_high = ceil(((seqstart - my_orf_coords) + 4) / 3)
                    
                    for i in range(1, 4):
                        Mut = MutableSeq(record.seq)
                        Mut[CutSiteCoords + (i - 1):seqstart + (i + 1)] = stoptag
                        Mut_AA = Mut[my_orf_coords:].translate(table=12)
                        
                        if Mut_AA[AA_low:AA_high].find('**') != -1:
                            GRIT = Seq(BsmBI_Up + Mut[CutSiteCoords + (i - 1) - HA_length:seqstart + (i - 1) + HA_length] + sequence + BsmBI_Down)
                            GRIT_List.append(GRIT)
                            
                            mid = len(GRIT) // 2
                            FW_Primer_List.append(Seq(GRIT[:mid + primer_len // 2]))
                            RV_Primer_List.append(Seq(GRIT[mid - primer_len // 2:]).reverse_complement())
                            Tm_List.append(mt.Tm_Wallace(GRIT[mid - primer_len // 2:mid + primer_len // 2]))
                            BindingRegion.append(GRIT[mid - primer_len // 2:mid + primer_len // 2])
                            
                            found = True
                            found_in_any_record = True
                            break
                
                if found:
                    break # Break out of FASTA loop if we found the frame
        
        if not found_in_any_record:
            print(f"Warning: gRNA '{sequence}' not found or could not generate valid in-frame mutation.")
            GRIT_List.append(np.nan)
            FW_Primer_List.append(np.nan)
            RV_Primer_List.append(np.nan)
            Tm_List.append(np.nan)
            BindingRegion.append(np.nan)

    # Write lists to dataframe 
    df['GRIT_Fragment'] = GRIT_List
    df['Forward_Primer'] = FW_Primer_List
    df['Reverse_Primer'] = RV_Primer_List
    df['Melting_Temp'] = Tm_List
    df['Binding_Region'] = BindingRegion
    
    # Output dataframe as csv file
    print(f"Writing results to: {args.output_csv}")
    df.to_csv(args.output_csv, index=False)
    print("Done!")

if __name__ == "__main__":
    main()
