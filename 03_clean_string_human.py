import pandas as pd
import os

# ===============================
TAX_ID = "9606"
CONF_THRESHOLD = 900
# ===============================

input_path = f"../data/raw/string/human/{TAX_ID}.protein.links.v12.0.txt"
output_dir = "../data/interim/cleaned_string/"
os.makedirs(output_dir, exist_ok=True)

df = pd.read_csv(input_path, sep=" ")

print("Original interactions:", len(df))

# Filter high confidence interactions
df = df[df["combined_score"] >= CONF_THRESHOLD].copy()

print("After confidence filter:", len(df))

output_path = output_dir + "human_string_clean.csv"
df.to_csv(output_path, index=False)

print("Saved:", output_path)
