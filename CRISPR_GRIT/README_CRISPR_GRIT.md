# CRISPR-GRIT Guide Generator

This Python script automates the design of CRISPR-GRIT repair templates. It aligns a list of guide RNAs (gRNAs) to a target gene FASTA file, calculates the precise cut site, inserts an in-frame mutation (e.g., stop codons), and extracts upstream and downstream homology arms. The final output is a complete repair fragment flanked by user-defined Golden Gate cloning overhangs (e.g., BsmBI), along with automatically generated forward and reverse primers for PCR amplification.

### Prerequisites
* Python 3.x
* Required Python packages: `biopython`, `pandas`, `numpy`

You can install the required dependencies using pip:
```bash
pip install biopython pandas numpy
```

### Inputs
1. **gRNA CSV File:** A comma-separated file containing your guide RNAs. It MUST contain a column named `Sequence` with the 20 bp gRNA sequences.
2. **Target Gene FASTA:** A standard `.fasta` file containing the genomic sequence of your target gene.

### Outputs
* **Output CSV File:** An augmented version of your input CSV that includes the following new columns:
  * `CRISPR-GRIT_Fragment`: The full sequence of the repair template.
  * `Forward_Primer` & `Reverse_Primer`: Primers designed to amplify the fragment.
  * `Melting_Temp`: Calculated Tm for the primer binding region.
  * `Binding_Region`: The specific overlap region used for primer design.

### Usage
This script is executed via the command line. It requires several arguments to define the homology arms, cloning overhangs, and primer overlap lengths.

**Basic Command Structure:**
```bash
python generate_grit_guides.py -i <input.csv> -f <gene.fasta> -o <output.csv> \
  --up_overhang <5bp_seq> --down_overhang <5bp_seq> \
  --ha_length <int> --overlap_length <int>
```

**Required Arguments:**
* `-i`, `--input_csv`: Path to the input gRNA CSV file.
* `-f`, `--fasta`: Path to the target gene FASTA file.
* `-o`, `--output_csv`: Path for the resulting output CSV file.
* `--up_overhang`: 5 bp upstream Golden Gate overhang (e.g., `GTCCA`).
* `--down_overhang`: 5 bp downstream Golden Gate overhang (e.g., `GTTTT`).
* `--ha_length`: Desired length of the homology arms in base pairs (e.g., `50`).
* `--overlap_length`: Desired oligo overlap length for primer design (e.g., `20`).

**Optional Arguments:**
* `--orf_start`: The nucleotide coordinate where the Open Reading Frame (ORF) starts in your FASTA file. Used to ensure stop codons are inserted in-frame. (Default: `100`)
* `--insert_seq`: The specific sequence to insert at the cut site. (Default: `taataa`)

### Example
```bash
python generate_grit_guides.py -i my_guides.csv -f ADE2.fasta -o GRIT_output.csv --up_overhang GTCCA --down_overhang GTTTT --ha_length 50 --overlap_length 20
```
