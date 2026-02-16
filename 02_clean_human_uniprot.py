from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import pandas as pd

# -----------------------------------
# STEP 1: Load FASTA
# -----------------------------------

fasta_path = "../data/raw/uniprot/human/human_uniprot_reviewed.fasta"

records = []

for record in SeqIO.parse(fasta_path, "fasta"):
    header = record.description
    sequence = str(record.seq)

    # Extract UniProt ID
    uniprot_id = header.split("|")[1] if "|" in header else header.split()[0]

    records.append({
        "uniprot_id": uniprot_id,
        "header": header,
        "sequence": sequence,
        "length": len(sequence)
    })

df_human = pd.DataFrame(records)

print("Total sequences:", len(df_human))


# -----------------------------------
# STEP 2: Remove invalid amino acids
# -----------------------------------

valid_amino_acids = set("ARNDCQEGHILKMFPSTWYV")

def is_valid_sequence(seq):
    return set(seq).issubset(valid_amino_acids)

df_human["is_valid"] = df_human["sequence"].apply(is_valid_sequence)

print("Invalid sequences count:", (~df_human["is_valid"]).sum())

df_human = df_human[df_human["is_valid"] == True].copy()

print("Remaining sequences after removing invalid ones:", len(df_human))


# -----------------------------------
# STEP 3: Remove duplicate sequences
# -----------------------------------

before_duplicates = len(df_human)

df_human = df_human.drop_duplicates(subset=["sequence"]).copy()

after_duplicates = len(df_human)

print("Duplicates removed:", before_duplicates - after_duplicates)
print("Remaining sequences after duplicate removal:", after_duplicates)


# -----------------------------------
# STEP 4: Length distribution
# -----------------------------------

print("Minimum length:", df_human["length"].min())
print("Maximum length:", df_human["length"].max())
print("Average length:", df_human["length"].mean())

print("\nLength statistics:")
print(df_human["length"].describe())


# -----------------------------------
# STEP 5: Length filtering
# -----------------------------------

before_length_filter = len(df_human)

df_human = df_human[
    (df_human["length"] >= 50) &
    (df_human["length"] <= 1500)
].copy()

after_length_filter = len(df_human)

print("Removed due to length filtering:", before_length_filter - after_length_filter)
print("Final clean sequence count:", after_length_filter)


# -----------------------------------
# STEP 6: Save cleaned dataset
# -----------------------------------

output_fasta = "../data/interim/cleaned_sequences/human_clean.fasta"
output_csv = "../data/interim/cleaned_sequences/human_clean.csv"

records_to_save = []

for _, row in df_human.iterrows():
    record = SeqRecord(
        Seq(row["sequence"]),
        id=row["uniprot_id"],
        description=""
    )
    records_to_save.append(record)

SeqIO.write(records_to_save, output_fasta, "fasta")
df_human.to_csv(output_csv, index=False)

print("Cleaned FASTA saved successfully!")
print("Cleaned CSV saved successfully!")
