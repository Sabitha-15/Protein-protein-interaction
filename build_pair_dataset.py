import os
import numpy as np
import pandas as pd

# =====================================================
# ROOT DIRECTORY (AUTO-DETECT PROJECT ROOT)
# =====================================================

BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

FEATURE_PATH = os.path.join(BASE_DIR, "data", "processed", "protein_features.npy")
DATASET_FOLDER = os.path.join(BASE_DIR, "data", "processed", "final_dataset")
OUTPUT_FOLDER = os.path.join(BASE_DIR, "data", "processed", "ml_ready")

os.makedirs(OUTPUT_FOLDER, exist_ok=True)

# =====================================================
# LOAD PROTEIN FEATURES
# =====================================================

print("🔹 Loading protein features...")

if not os.path.exists(FEATURE_PATH):
    print(f"❌ Protein features file not found at: {FEATURE_PATH}")
    exit()

raw_features = np.load(FEATURE_PATH, allow_pickle=True).item()

# Clean UniProt IDs (sp|P07001|PNTA_ECOLI → P07001)
protein_features = {}

for key, value in raw_features.items():
    if "|" in key:
        clean_id = key.split("|")[1]
    else:
        clean_id = key
    protein_features[clean_id] = value

print("✅ Total proteins loaded:", len(protein_features))

# =====================================================
# FUNCTION TO BUILD DATASET
# =====================================================

def build_dataset(excel_filename, species_name):

    file_path = os.path.join(DATASET_FOLDER, excel_filename)

    if not os.path.exists(file_path):
        print(f"❌ Dataset file not found: {file_path}")
        return  # Skip processing

    print(f"\n🔹 Processing {species_name} dataset...")
    print("Reading file:", file_path)

    df = pd.read_excel(file_path)

    required_cols = ["protein1_uniprot", "protein2_uniprot", "label"]
    for col in required_cols:
        if col not in df.columns:
            print(f"❌ Missing required column '{col}' in {excel_filename}")
            return

    X = []
    y = []

    missing_count = 0
    missing_proteins = []

    for _, row in df.iterrows():

        p1 = str(row["protein1_uniprot"]).strip()
        p2 = str(row["protein2_uniprot"]).strip()
        label = int(row["label"])

        if p1 in protein_features and p2 in protein_features:

            f1 = protein_features[p1]
            f2 = protein_features[p2]

            # Simple concatenation (protein1 + protein2)
            pair_feature = np.concatenate([f1, f2])

            X.append(pair_feature)
            y.append(label)

        else:
            missing_count += 1
            missing_proteins.append((p1, p2))

    if len(X) == 0:
        print(f"❌ No valid protein pairs found for {species_name}. Skipping save.")
        return

    X = np.array(X)
    y = np.array(y)

    print("✅ Total valid pairs:", len(X))
    print("⚠ Missing proteins skipped:", missing_count)
    if missing_count > 0:
        print("⚠ Example missing pairs:", missing_proteins[:5])
    print("📊 Feature size per pair:", X.shape[1])

    # Save files
    np.save(os.path.join(OUTPUT_FOLDER, f"X_{species_name}.npy"), X)
    np.save(os.path.join(OUTPUT_FOLDER, f"y_{species_name}.npy"), y)

    print(f"💾 {species_name} dataset saved successfully.")


# =====================================================
# BUILD ECOli & HUMAN DATASETS
# =====================================================

build_dataset("ecoli_ppi_balanced_dataset.xlsx", "ecoli")
build_dataset("human_ppi_balanced_dataset.xlsx", "human")

print("\n🎉 ALL DATASETS READY FOR MACHINE LEARNING")
