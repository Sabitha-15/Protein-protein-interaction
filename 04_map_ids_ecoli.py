import pandas as pd
import os

# ==========================================
# STEP 1: Load cleaned STRING interactions
# ==========================================

string_path = "../data/interim/cleaned_sequences/ecoli_string_clean.csv"
df_string = pd.read_csv(string_path)

print("Loaded STRING high-confidence interactions:", len(df_string))


# ==========================================
# STEP 2: Load STRING alias file
# ==========================================

alias_path = "../data/raw/string/ecoli/511145.protein.aliases.v12.0.txt"
df_alias = pd.read_csv(alias_path, sep="\t")

# Keep only UniProt mappings
df_uniprot_alias = df_alias[df_alias["source"] == "UniProt_AC"].copy()

print("Total UniProt alias entries:", len(df_uniprot_alias))


# ==========================================
# STEP 3: Load cleaned sequences
# ==========================================

seq_path = "../data/interim/cleaned_sequences/ecoli_clean.csv"
df_sequences = pd.read_csv(seq_path)

clean_uniprot_ids = set(df_sequences["uniprot_id"])

# Keep only mappings present in cleaned sequences
df_uniprot_alias_filtered = df_uniprot_alias[
    df_uniprot_alias["alias"].isin(clean_uniprot_ids)
].copy()

print("Valid STRING → UniProt mappings:", len(df_uniprot_alias_filtered))


# ==========================================
# STEP 4: Build mapping dictionary
# ==========================================

mapping_dict = dict(
    zip(
        df_uniprot_alias_filtered["#string_protein_id"],
        df_uniprot_alias_filtered["alias"]
    )
)

print("Total unique STRING IDs mapped:", len(mapping_dict))


# ==========================================
# STEP 5: Map STRING IDs to UniProt
# ==========================================

df_string["protein1_uniprot"] = df_string["protein1"].map(mapping_dict)
df_string["protein2_uniprot"] = df_string["protein2"].map(mapping_dict)

before_mapping = len(df_string)

df_string = df_string.dropna(
    subset=["protein1_uniprot", "protein2_uniprot"]
).copy()

after_mapping = len(df_string)

print("Removed due to missing mapping:", before_mapping - after_mapping)
print("Mapped interaction pairs:", after_mapping)


# ==========================================
# STEP 6: Attach sequences
# ==========================================

sequence_dict = dict(
    zip(df_sequences["uniprot_id"], df_sequences["sequence"])
)

df_string["sequence1"] = df_string["protein1_uniprot"].map(sequence_dict)
df_string["sequence2"] = df_string["protein2_uniprot"].map(sequence_dict)

before_seq = len(df_string)

df_string = df_string.dropna(
    subset=["sequence1", "sequence2"]
).copy()

after_seq = len(df_string)

print("Removed due to missing sequences:", before_seq - after_seq)
print("Final positive interaction dataset size:", after_seq)


# ==========================================
# STEP 7: Add label column
# ==========================================

df_string["label"] = 1

print("Label column added.")


# ==========================================
# STEP 8: Save positive dataset
# ==========================================

output_dir = "../data/processed/interaction_pairs/"
os.makedirs(output_dir, exist_ok=True)

output_path = output_dir + "ecoli_positive_pairs.csv"

df_string.to_csv(output_path, index=False)

print("\n✅ Positive interaction dataset saved successfully!")
print("Saved at:", output_path)
