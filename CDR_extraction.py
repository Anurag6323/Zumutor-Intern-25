from Bio import SeqIO
import csv

# Define CDR ranges (0-based index, end-exclusive)
CDR1_RANGE = (30, 40)
CDR2_RANGE = (55, 65)
CDR3_RANGE = (95, 105)


# File Paths
input_fasta = "trem1_variants.fasta"    # Input FASTA
output_csv = "trem1_cdrs.csv"           # Output in .csv file

# Extraction
with open(output_csv, "w", newline='') as csvfile:
    writer = csv.writer(csvfile)
    writer.writerow(["Sequence_ID", "CDR1", "CDR2", "CDR3"])  # Header

    for record in SeqIO.parse(input_fasta, "fasta"):
        seq_id = record.id
        sequence = str(record.seq)

        # Check if sequence is long enough
        if len(sequence) < max(CDR1_RANGE[1], CDR2_RANGE[1], CDR3_RANGE[1]):
            print(f"Warning: {seq_id} is shorter than expected. Skipping.")
            continue

        # Extract CDRs
        cdr1 = sequence[CDR1_RANGE[0]:CDR1_RANGE[1]]
        cdr2 = sequence[CDR2_RANGE[0]:CDR2_RANGE[1]]
        cdr3 = sequence[CDR3_RANGE[0]:CDR3_RANGE[1]]

        writer.writerow([seq_id, cdr1, cdr2, cdr3])

print(f"Extraction complete! CDR data saved to '{output_csv}'")
