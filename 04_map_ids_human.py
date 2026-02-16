import pandas as pd
import os

TAX_ID = "9606"

# Paths
string_path = "../data/interim/cleaned_string/human_string_clean.csv"
alias_path = f"../data/raw/string/human/{TAX_ID}.protein.aliases.v12.0.txt"
seq_path = "../data/interim/cleaned_sequences/human_clean.csv"

df_string = pd.read_csv(string_path)
df_alias = pd.read_csv(alias_path, sep="\t")
df_sequences = pd.read_csv(seq_path)

# Keep UniProt mappings only
df_alias = df_alias[df_alias["source"] == "UniProt_AC"]

# Keep only UniProt IDs present in cleaned sequences
clean_ids = set(df_sequences["uniprot_id"])
df_alias = df_alias[df_alias["alias"].isin(clean_ids)]

# Create mapping
mapping_dict = dict(zip(
    df_alias["#string_protein_id"],
    df_alias["alias"]
))

# Map IDs
df_string["protein1_uniprot"] = df_string["protein1"].map(mapping_dict)
df_string["protein2_uniprot"] = df_string["protein2"].map(mapping_dict)

df_string = df_string.dropna(subset=["protein1_uniprot", "protein2_uniprot"])

# Attach sequences
seq_dict = dict(zip(df_sequences["uniprot_id"], df_sequences["sequence"]))

df_string["sequence1"] = df_string["protein1_uniprot"].map(seq_dict)
df_string["sequence2"] = df_string["protein2_uniprot"].map(seq_dict)

df_string = df_string.dropna(subset=["sequence1", "sequence2"])

df_string["label"] = 1

# Save
output_dir = "../data/processed/interaction_pairs/"
os.makedirs(output_dir, exist_ok=True)

output_path = output_dir + "human_positive_pairs.csv"
df_string.to_csv(output_path, index=False)

print("Saved:", output_path)
print("Total positive pairs:", len(df_string))
