from Bio import SeqIO
import csv
from collections import defaultdict

# === EDIT CDR RANGES HERE (1-based, inclusive) ===
CDR_L1_RANGE = (24, 38)   # Light chain CDR1
CDR_H2_RANGE = (50, 68)   # Heavy chain CDR2
CDR_H3_RANGE = (101, 110) # Heavy chain CDR3
# ==================================================

# Convert to 0-based for Python slicing (end-exclusive)
cdr_l1_range = (CDR_L1_RANGE[0] - 1, CDR_L1_RANGE[1])
cdr_h2_range = (CDR_H2_RANGE[0] - 1, CDR_H2_RANGE[1])
cdr_h3_range = (CDR_H3_RANGE[0] - 1, CDR_H3_RANGE[1])

# Input/Output files
input_fasta = "trem1_variants.fasta"
output_csv = "trem1_cdr_pairs.csv"

# Dictionary to group light and heavy chains by mutation ID
mutations = defaultdict(dict)

# Load sequences and separate by chain
for record in SeqIO.parse(input_fasta, "fasta"):
    seq_id = record.id
    sequence = str(record.seq)

    # Expected format: mutX_L or mutX_H
    if seq_id.endswith("_L"):
        mut_id = seq_id[:-2]
        mutations[mut_id]["L"] = sequence
    elif seq_id.endswith("_H"):
        mut_id = seq_id[:-2]
        mutations[mut_id]["H"] = sequence
    else:
        print(f"Warning: Unrecognized sequence ID format: {seq_id}")

# Write CDR data
with open(output_csv, "w", newline='') as csvfile:
    writer = csv.writer(csvfile)
    writer.writerow(["Mutation_ID", "CDR_L1", "CDR_H2", "CDR_H3"])  # Header

    for mut_id, chains in mutations.items():
        light_chain = chains.get("L")
        heavy_chain = chains.get("H")

        if not light_chain or not heavy_chain:
            print(f"Warning: Missing chain for {mut_id}, skipping.")
            continue

        try:
            cdr_l1 = light_chain[cdr_l1_range[0]:cdr_l1_range[1]]
            cdr_h2 = heavy_chain[cdr_h2_range[0]:cdr_h2_range[1]]
            cdr_h3 = heavy_chain[cdr_h3_range[0]:cdr_h3_range[1]]
            writer.writerow([mut_id, cdr_l1, cdr_h2, cdr_h3])
        except IndexError:
            print(f"Warning: Sequence too short in {mut_id}, skipping.")

print(f"\n✅ Done! Extracted CDRs saved to '{output_csv}'")
