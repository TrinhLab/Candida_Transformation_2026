# gRNA Library FASTQ Counting

This script takes an input sequencing FASTQ file and a gRNA library CSV file and outputs gRNA read counts. 

### Prerequisites
* Python 3
* Pandas (`pip install pandas`)

### Usage
Run the script from the terminal by providing the paths to your FASTQ file, your metadata CSV, and your desired output CSV name:

python gRNALibraryFastQCounting.py -i <input.fastq> -m <metadata.csv> -o <output.csv>

**Example:**
python gRNALibraryFastQCounting.py -i CaVipLibPool_Rep3Joined.fastq -m CaVipLibraryMetadata.csv -o Rep3_Counts.csv
