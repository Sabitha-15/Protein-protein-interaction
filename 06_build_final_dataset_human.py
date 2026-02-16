import pandas as pd
import os

pos_path = "../data/processed/interaction_pairs/human_positive_pairs.csv"
neg_path = "../data/processed/interaction_pairs/human_negative_pairs.csv"

df_pos = pd.read_csv(pos_path)
df_neg = pd.read_csv(neg_path)

df_final = pd.concat([df_pos, df_neg], ignore_index=True)

df_final = df_final.sample(frac=1, random_state=42).reset_index(drop=True)

output_dir = "../data/processed/final_dataset/"
os.makedirs(output_dir, exist_ok=True)

output_path = output_dir + "human_ppi_balanced_dataset.csv"
df_final.to_csv(output_path, index=False)

print("Saved:", output_path)
print("Total samples:", len(df_final))
