import pandas as pd
import random
import os

# -----------------------------------
# STEP 1: Load positive dataset
# -----------------------------------

df_positive = pd.read_csv("../data/processed/interaction_pairs/ecoli_positive_pairs.csv")

print("Positive pairs:", len(df_positive))

# -----------------------------------
# STEP 2: Get unique protein list
# -----------------------------------

proteins = list(set(df_positive["protein1_uniprot"]).union(
                set(df_positive["protein2_uniprot"])))

print("Total unique proteins:", len(proteins))

# -----------------------------------
# STEP 3: Create set of positive pairs
# (order-independent)
# -----------------------------------

positive_pairs = set()

for _, row in df_positive.iterrows():
    pair = tuple(sorted([row["protein1_uniprot"], row["protein2_uniprot"]]))
    positive_pairs.add(pair)

print("Total unique positive pairs:", len(positive_pairs))

# -----------------------------------
# STEP 4: Generate negative pairs
# -----------------------------------

negative_pairs = set()

target_negatives = len(positive_pairs)

while len(negative_pairs) < target_negatives:
    p1, p2 = random.sample(proteins, 2)
    pair = tuple(sorted([p1, p2]))

    if pair not in positive_pairs:
        negative_pairs.add(pair)

print("Generated negative pairs:", len(negative_pairs))

# -----------------------------------
# STEP 5: Convert to DataFrame
# -----------------------------------

df_negative = pd.DataFrame(list(negative_pairs),
                           columns=["protein1_uniprot", "protein2_uniprot"])

df_negative["label"] = 0

# -----------------------------------
# STEP 6: Attach sequences
# -----------------------------------

df_sequences = pd.read_csv("../data/interim/cleaned_sequences/ecoli_clean.csv")

sequence_dict = dict(zip(df_sequences["uniprot_id"], df_sequences["sequence"]))

df_negative["sequence1"] = df_negative["protein1_uniprot"].map(sequence_dict)
df_negative["sequence2"] = df_negative["protein2_uniprot"].map(sequence_dict)

print("Negative dataset ready.")

# -----------------------------------
# STEP 7: Save negative dataset
# -----------------------------------

output_path = "../data/processed/interaction_pairs/ecoli_negative_pairs.csv"
os.makedirs("../data/processed/interaction_pairs/", exist_ok=True)

df_negative.to_csv(output_path, index=False)

print("Negative pairs saved successfully!")
