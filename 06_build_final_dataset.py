import pandas as pd
import os

# ==========================================
# STEP 1: Load positive and negative datasets
# ==========================================

pos_path = "../data/processed/interaction_pairs/ecoli_positive_pairs.csv"
neg_path = "../data/processed/interaction_pairs/ecoli_negative_pairs.csv"

df_pos = pd.read_csv(pos_path)
df_neg = pd.read_csv(neg_path)

print("Positive samples:", len(df_pos))
print("Negative samples:", len(df_neg))


# ==========================================
# STEP 2: Keep only necessary columns
# ==========================================

columns_needed = [
    "protein1_uniprot",
    "protein2_uniprot",
    "sequence1",
    "sequence2",
    "label"
]

df_pos = df_pos[columns_needed]
df_neg = df_neg[columns_needed]


# ==========================================
# STEP 3: Combine datasets
# ==========================================

df_final = pd.concat([df_pos, df_neg], ignore_index=True)

print("Total combined dataset:", len(df_final))


# ==========================================
# STEP 4: Shuffle dataset
# ==========================================

df_final = df_final.sample(frac=1, random_state=42).reset_index(drop=True)

print("Dataset shuffled.")


# ==========================================
# STEP 5: Save final dataset
# ==========================================

output_dir = "../data/processed/final_dataset/"
os.makedirs(output_dir, exist_ok=True)

output_path = output_dir + "ecoli_ppi_balanced_dataset.csv"

df_final.to_csv(output_path, index=False)

print("\n✅ Final ML-ready dataset saved successfully!")
print("Saved at:", output_path)
