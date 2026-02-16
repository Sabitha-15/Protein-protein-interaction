import pandas as pd
import random
import os

pos_path = "../data/processed/interaction_pairs/human_positive_pairs.csv"
df_pos = pd.read_csv(pos_path)

print("Positive pairs:", len(df_pos))

# Unique proteins
proteins = set(df_pos["protein1_uniprot"]).union(
    set(df_pos["protein2_uniprot"])
)

proteins = list(proteins)

# Unique positive pairs
pos_pairs = set(
    tuple(sorted([a, b]))
    for a, b in zip(df_pos["protein1_uniprot"], df_pos["protein2_uniprot"])
)

num_neg = len(pos_pairs)

neg_pairs = set()

while len(neg_pairs) < num_neg:
    p1, p2 = random.sample(proteins, 2)
    pair = tuple(sorted([p1, p2]))
    if pair not in pos_pairs:
        neg_pairs.add(pair)

print("Generated negative pairs:", len(neg_pairs))

df_neg = pd.DataFrame(list(neg_pairs),
                      columns=["protein1_uniprot", "protein2_uniprot"])

# Attach sequences
seq_dict = dict(zip(
    df_pos["protein1_uniprot"],
    df_pos["sequence1"]
))

df_neg["sequence1"] = df_neg["protein1_uniprot"].map(seq_dict)
df_neg["sequence2"] = df_neg["protein2_uniprot"].map(seq_dict)

df_neg["label"] = 0

output_dir = "../data/processed/interaction_pairs/"
os.makedirs(output_dir, exist_ok=True)

output_path = output_dir + "human_negative_pairs.csv"
df_neg.to_csv(output_path, index=False)

print("Saved:", output_path)
