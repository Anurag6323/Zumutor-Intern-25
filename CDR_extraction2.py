from Bio import SeqIO
import re
from collections import defaultdict
import pandas as pd

# === Define CDR ranges (1-based, inclusive) ===
CDR_L1_RANGE = (24, 38)
CDR_H2_RANGE = (50, 68)
CDR_H3_RANGE = (101, 110)
# ==============================================

# Convert to 0-based Python slicing
cdr_l1_range = (CDR_L1_RANGE[0] - 1, CDR_L1_RANGE[1])
cdr_h2_range = (CDR_H2_RANGE[0] - 1, CDR_H2_RANGE[1])
cdr_h3_range = (CDR_H3_RANGE[0] - 1, CDR_H3_RANGE[1])

input_fasta = "trem1_variants.fasta"
output_excel = "trem1_cdrs_with_mutations.xlsx"

mutations = defaultdict(dict)

# Load and group sequences
for record in SeqIO.parse(input_fasta, "fasta"):
    header = record.description
    sequence = str(record.seq)

    match = re.match(r"(mut\d+)_([LH])\|([^|]*)\|", header)
    if not match:
        print(f"Warning: Unrecognized header format: {header}")
        continue

    mut_id, chain_type, mut_list_str = match.groups()
    mut_list = [m.strip() for m in mut_list_str.split(",") if m.strip()]
    mutations[mut_id][chain_type] = {
        "seq": sequence,
        "mutations": mut_list
    }

def extract_and_highlight(seq, mut_list, cdr_range):
    cdr_seq = list(seq[cdr_range[0]:cdr_range[1]])
    cdr_muts = []

    for mut in mut_list:
        match = re.match(r"([A-Z])(\d+)([A-Z])", mut)
        if not match:
            continue
        wt, pos_str, new = match.groups()
        pos = int(pos_str) - 1
        if cdr_range[0] <= pos < cdr_range[1]:
            rel_pos = pos - cdr_range[0]
            cdr_seq[rel_pos] = f"[{cdr_seq[rel_pos]}]"
            cdr_muts.append(mut)

    return ''.join(cdr_seq), ', '.join(cdr_muts)

# Build DataFrame rows
rows = []

for mut_id, chains in mutations.items():
    if "L" not in chains or "H" not in chains:
        print(f"Warning: Missing L/H chain for {mut_id}, skipping.")
        continue

    light = chains["L"]
    heavy = chains["H"]

    try:
        cdr_l1, muts_l1 = extract_and_highlight(light["seq"], light["mutations"], cdr_l1_range)
        cdr_h2, muts_h2 = extract_and_highlight(heavy["seq"], heavy["mutations"], cdr_h2_range)
        cdr_h3, muts_h3 = extract_and_highlight(heavy["seq"], heavy["mutations"], cdr_h3_range)

        rows.append({
            "Mutation_ID": mut_id,
            "CDR_L1 (light)": cdr_l1,
            "CDR_H2 (heavy)": cdr_h2,
            "CDR_H3 (heavy)": cdr_h3,
            "CDR_L1_Muts": muts_l1,
            "CDR_H2_Muts": muts_h2,
            "CDR_H3_Muts": muts_h3
        })
    except IndexError:
        print(f"Warning: Sequence too short in {mut_id}, skipping.")

# Export to Excel
df = pd.DataFrame(rows)
df.to_excel(output_excel, index=False)

print(f"\n✅ Done! Excel file saved to '{output_excel}'")
