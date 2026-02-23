import os
import numpy as np
from Bio import SeqIO
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from collections import Counter
import itertools
import glob

# =====================================================
# ROOT DIRECTORY
# =====================================================

BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

# Look for all .fasta, .fa, .faa files in uniprot folder and subfolders
DATA_PATH_PATTERN = os.path.join(BASE_DIR, "data", "raw", "uniprot", "**", "*.*a*")
fasta_files = glob.glob(DATA_PATH_PATTERN, recursive=True)

OUTPUT_PATH = os.path.join(BASE_DIR, "data", "processed", "protein_features.npy")

# =====================================================
# AMINO ACIDS
# =====================================================

valid_aa = set("ACDEFGHIKLMNPQRSTVWY")
amino_acids = "ACDEFGHIKLMNPQRSTVWY"
dipeptides = [''.join(p) for p in itertools.product(amino_acids, repeat=2)]

# =====================================================
# FEATURE FUNCTIONS
# =====================================================

def compute_aac(seq):
    length = len(seq)
    count = Counter(seq)
    return np.array([count[aa] / length for aa in amino_acids])

def compute_dipeptide(seq):
    total = len(seq) - 1
    freq = dict.fromkeys(dipeptides, 0)

    for i in range(total):
        pair = seq[i:i+2]
        if pair in freq:
            freq[pair] += 1

    return np.array([freq[p] / total for p in dipeptides])

def compute_physico(seq):
    analysed = ProteinAnalysis(seq)
    helix, turn, sheet = analysed.secondary_structure_fraction()

    return np.array([
        analysed.molecular_weight(),
        analysed.aromaticity(),
        analysed.instability_index(),
        analysed.gravy(),
        analysed.isoelectric_point(),
        helix,
        turn,
        sheet
    ])

# =====================================================
# MAIN
# =====================================================

protein_features = {}

# Debug print to check which files are found
print("Looking for FASTA files under:", os.path.join(BASE_DIR, "data", "raw", "uniprot"))
print("Found files:", fasta_files)

if len(fasta_files) == 0:
    print("❌ No FASTA files found.")
    exit()

for fasta_file in fasta_files:
    print("\n📂 Reading:", fasta_file)

    for record in SeqIO.parse(fasta_file, "fasta"):
        raw_seq = str(record.seq)

        # Clean non-standard amino acids
        seq = "".join([aa for aa in raw_seq if aa in valid_aa])

        if len(seq) < 2:
            continue

        pid = record.id

        aac = compute_aac(seq)
        dipep = compute_dipeptide(seq)
        phys = compute_physico(seq)

        feature_vector = np.concatenate([aac, dipep, phys])
        protein_features[pid] = feature_vector

# Create processed folder if it doesn't exist
os.makedirs(os.path.dirname(OUTPUT_PATH), exist_ok=True)

# Save features
np.save(OUTPUT_PATH, protein_features)

print("\n✅ DONE")
print("Total proteins:", len(protein_features))
print("Feature size per protein:", len(feature_vector))
print("Saved to:", OUTPUT_PATH)
