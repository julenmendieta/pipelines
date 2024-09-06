import gzip
from collections import defaultdict
from collections import Counter

###
# Script to check in a Perturb-seq fastq file the presence of 10X barcodes 
# (cell IDs) and CRISPR gRNAs 
 
def process_fastq_files(read1_file, read2_file):
    """
    Processes paired-end FASTQ files for Pertub-seq data.

    Args:
        read1_file (str): Path to the Read1 FASTQ file.
        read2_file (str): Path to the Read2 FASTQ file.

    Returns:
        tuple: A tuple containing two dictionaries:
            - cell_barcode_counts: A dictionary mapping cell barcodes to counts.
            - gRNA_counts: A dictionary mapping gRNA sequences to counts.
    """

    cell_barcode_counts = defaultdict(int)
    gRNA_counts = defaultdict(int)
    cell_barcode_to_gRNA = defaultdict(list)
    gRNA_to_cell_barcode = defaultdict(list)

    with gzip.open(read1_file, 'rt') as f1, gzip.open(read2_file, 'rt') as f2:
        while True:
            # Read four lines from each file
            r1_header = f1.readline().strip()
            r1_seq = f1.readline().strip()[0:16]  # Extract cell barcode
            f1.readline()
            f1.readline()

            r2_header = f2.readline().strip()
            r2_seq = f2.readline().strip()[31:51]  # Extract gRNA
            f2.readline()
            f2.readline()

            if not r1_header or not r2_header:
                break

            cell_barcode_counts[r1_seq] += 1
            gRNA_counts[r2_seq] += 1
            cell_barcode_to_gRNA[r1_seq] += [r2_seq]
            gRNA_to_cell_barcode[r2_seq] += [r1_seq]

    return cell_barcode_counts, gRNA_counts, cell_barcode_to_gRNA, gRNA_to_cell_barcode

# Example usage


read1_file = "Perturb-Seq_LSK_d7_CRISPR_S5_R1_001.fastq.gz"
read2_file = "Perturb-Seq_LSK_d7_CRISPR_S5_R2_001.fastq.gz"


givenGuides = ["ATGTTGCAGTTCGGCTCGAT", "ACGTGTAAGGCGAACGCCTT", "GACTCCGGGTACTAAATGTC",
               "CCGCGCCGTTAGGGAACGAG", "AGGACGGTAGAAGTTTCAGG", "GCATCCATCTTCATTCACAG",
               "GGATCATCAAGACTCCCCGG", "AGAAAGGGCGGCGATCAAGG", "GCTGGTGCGCGAGTTCAGCG",
               "GACAGCGTTCACATTCCAGG", "GGGAGAGGTCCACAGCGCCG", "TCCTCACTCTTGGATAACAG",
               "CAATGCCAGGCACTACAATG", "CTGTCACTTTAATAGCTCTG", "GAATGCGCCCTAAATCACTG",
               "CATGGACGATAACGACCACA", "GTTCGTTGTACAGAATGTGA", "ACAGCAGCTCTATCGCCACA",
               "TGATAAGGGAAGCACATCCG", "CTTGATGGAACGGTTCAATG", "TCATCCAGTCGAGAAACCGG",
               "CCAGCGGGTGAAATACACCA", "CGCCGCAGCGAATAATTCGG", "TCCGGAAGCGAAGTTTAAAG",
               "TCTGCATCATCATCAAACTG"]

id1 = read1_file.split('_R1_001.fastq.gz')[0]
print(id1)
cell_barcode_counts, gRNA_counts, cell_barcode_to_gRNA, gRNA_to_cell_barcode = process_fastq_files(read1_file, read2_file)

# Sort dictionary by value
cell_barcode_counts = {k: v for k, v in sorted(cell_barcode_counts.items(), 
                                               key=lambda item: item[1])[::-1]}
gRNA_counts = {k: v for k, v in sorted(gRNA_counts.items(), 
                                               key=lambda item: item[1])[::-1]}

# So top N barcodes and gRNAs
topN = 100
#for c in list(cell_barcode_counts.keys())[:topN]:
#    print(c, cell_barcode_counts[c])
for c in list(gRNA_counts.keys())[:topN]:
    print(c, gRNA_counts[c])
print()

#[c for c in set(cell_barcode_to_gRNA["TGCGACGCTCAGACAG"]) if c in givenGuides]



#a = Counter(gRNA_to_cell_barcode['GACTCCGGGTACTAAATGTC'])
#a = {k: v for k, v in sorted(a.items(), key=lambda item: item[1])[::-1]}
#for c in list(a.keys())[:topN]:
#    print(c, a[c])