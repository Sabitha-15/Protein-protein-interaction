from Bio import SeqIO
import pandas as pd

# Path to your ecoli fasta file
fasta_path = fasta_path = "../data/raw/uniprot/ecoli/ecoli_uniprot_reviewed.fasta"

records = []

for record in SeqIO.parse(fasta_path, "fasta"):
    header = record.description
    sequence = str(record.seq)
    
    # Extract UniProt ID (second field between | )
    uniprot_id = header.split("|")[1] if "|" in header else header.split()[0]
    
    records.append({
        "uniprot_id": uniprot_id,
        "header": header,
        "sequence": sequence,
        "length": len(sequence)
    })


df_ecoli = pd.DataFrame(records)

print("Total sequences:", len(df_ecoli))
df_ecoli.head()

# -----------------------------------
# STEP 2: Remove invalid amino acids
# -----------------------------------

# Standard 20 amino acids
valid_amino_acids = set("ARNDCQEGHILKMFPSTWYV")

def is_valid_sequence(seq):
    return set(seq).issubset(valid_amino_acids)

# Check invalid sequences
df_ecoli["is_valid"] = df_ecoli["sequence"].apply(is_valid_sequence)

print("Invalid sequences count:", (~df_ecoli["is_valid"]).sum())

# Keep only valid sequences
df_ecoli = df_ecoli[df_ecoli["is_valid"] == True].copy()

print("Remaining sequences after removing invalid ones:", len(df_ecoli))
# -----------------------------------
# STEP 3: Remove duplicate sequences
# -----------------------------------

before_duplicates = len(df_ecoli)

df_ecoli = df_ecoli.drop_duplicates(subset=["sequence"]).copy()

after_duplicates = len(df_ecoli)

print("Duplicates removed:", before_duplicates - after_duplicates)
print("Remaining sequences after duplicate removal:", after_duplicates)

# -----------------------------------
# STEP 4: Length distribution
# -----------------------------------

print("Minimum length:", df_ecoli["length"].min())
print("Maximum length:", df_ecoli["length"].max())
print("Average length:", df_ecoli["length"].mean())

print("\nLength statistics:")
print(df_ecoli["length"].describe())

# -----------------------------------
# STEP 5: Length filtering
# -----------------------------------

before_length_filter = len(df_ecoli)

df_ecoli = df_ecoli[
    (df_ecoli["length"] >= 50) & 
    (df_ecoli["length"] <= 1500)
].copy()

after_length_filter = len(df_ecoli)

print("Removed due to length filtering:", before_length_filter - after_length_filter)
print("Final clean sequence count:", after_length_filter)

# -----------------------------------
# STEP 6: Save cleaned dataset
# -----------------------------------

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# Save cleaned FASTA
output_fasta = "../data/interim/cleaned_sequences/ecoli_clean.fasta"

records_to_save = []

for _, row in df_ecoli.iterrows():
    record = SeqRecord(
        Seq(row["sequence"]),
        id=row["uniprot_id"],
        description=""
    )
    records_to_save.append(record)

SeqIO.write(records_to_save, output_fasta, "fasta")

print("Cleaned FASTA saved successfully!")

# Save CSV version
output_csv = "../data/interim/cleaned_sequences/ecoli_clean.csv"
df_ecoli.to_csv(output_csv, index=False)

print("Cleaned CSV saved successfully!")
