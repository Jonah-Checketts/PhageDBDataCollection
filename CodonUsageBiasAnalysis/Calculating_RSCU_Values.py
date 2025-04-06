#!/usr/bin/env python3

import os
import sys
from Bio import SeqIO
from collections import defaultdict
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# Codon usage bias table:
CODON_TABLE = {
    'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
    'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
    'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
    'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
    'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
    'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
    'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
    'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
    'TAT': 'Y', 'TAC': 'Y', 'CAT': 'H', 'CAC': 'H',
    'CAA': 'Q', 'CAG': 'Q', 'AAT': 'N', 'AAC': 'N',
    'AAA': 'K', 'AAG': 'K', 'GAT': 'D', 'GAC': 'D',
    'GAA': 'E', 'GAG': 'E', 'TGT': 'C', 'TGC': 'C',
    'TGG': 'W', 'CGT': 'R', 'CGC': 'R', 'CGA': 'R',
    'CGG': 'R', 'AGT': 'S', 'AGC': 'S', 'AGA': 'R',
    'AGG': 'R', 'GGT': 'G', 'GGC': 'G', 'GGA': 'G',
    'GGG': 'G'
}

def clean_cds(filepath):
    """
    For bacterial sequences which may not always be multiples of 3.
    Removes invalid CDS from the file in-place.
    """
    good, bad = [], []
    for rec in SeqIO.parse(filepath, "fasta"):
        if len(rec.seq) % 3 == 0:
            good.append(rec)
        else:
            bad.append(rec.id)
    SeqIO.write(good, filepath, "fasta")
    print(f"Removed {len(bad)} invalid CDS: {bad[:10]}")

def Get_Synonymous_Codons(codon):
    count = 0
    synonymous_aa = CODON_TABLE[codon]  # The amino acid for this codon
    for c in CODON_TABLE:
        if CODON_TABLE[c] == synonymous_aa:
            count += 1
    return count

def Calculate_RSCU(codon_counts, amino_acid_counts):
    RSCU_data = {}
    # RSCU = (Observed codon count) / ( (Total codons for that AA) / (Number of synonymous codons for that AA) )
    for cd in CODON_TABLE:
        if amino_acid_counts[CODON_TABLE[cd]] == 0:
            # Avoid division by zero if there's an amino acid never observed
            RSCU_data[cd] = 0
            continue
        R_val = codon_counts[cd] / (
            amino_acid_counts[CODON_TABLE[cd]] / Get_Synonymous_Codons(cd)
        )
        RSCU_data[cd] = R_val
    return RSCU_data

def Parse_CDS(folder_or_file):
    """
    If given a folder, parses all .fasta inside it.
    If given a single file, just parses that file.
    Collects codon usage data and returns an RSCU dictionary.
    """
    codon_counts = defaultdict(int)
    amino_acid_counts = defaultdict(int)

    if os.path.isdir(folder_or_file):
        # We'll parse every .fasta in the folder
        for fname in os.listdir(folder_or_file):
            if fname.endswith('.fasta'):
                file_path = os.path.join(folder_or_file, fname)
                #clean_cds(file_path)  # Optionally clean
                with open(file_path, "r") as f:
                    for record in SeqIO.parse(f, "fasta"):
                        accumulate_codon_counts(record, codon_counts, amino_acid_counts)
    else:
        # It's a single file
        #clean_cds(folder_or_file)  # Optionally clean
        with open(folder_or_file, "r") as f:
            for record in SeqIO.parse(f, "fasta"):
                accumulate_codon_counts(record, codon_counts, amino_acid_counts)

    return Calculate_RSCU(codon_counts, amino_acid_counts)

def accumulate_codon_counts(record, codon_counts, amino_acid_counts):
    """Helper to scan a sequence record in 3-nt windows."""
    seq_len = len(record.seq)
    if seq_len % 3 != 0:
        # You could skip or raise an error
        return
    for i in range(0, seq_len, 3):
        codon = str(record.seq[i:i+3]).upper()
        # Skip stop codons
        if codon in ("TGA","TAA","TAG"):
            continue
        codon_counts[codon] += 1
        amino_acid_counts[CODON_TABLE[codon]] += 1

def main():
    # -----------------------------
    # 1) Ask user for bacterial info
    # -----------------------------
    bacteria_path = input("Enter the path to your bacterial sequence folder/file: ").strip()
    bacteria_gc_str = input("Enter the GC% for the bacterial genome (e.g. 66.4): ").strip()
    try:
        bacteria_gc = float(bacteria_gc_str)
    except ValueError:
        bacteria_gc = None  # or default to something

    # -----------------------------
    # 2) Ask how many phage clusters
    # -----------------------------
    while True:
        try:
            num_clusters = int(input("How many phage clusters do you want to analyze? "))
            break
        except ValueError:
            print("Please enter an integer.")

    cluster_info = []
    for i in range(num_clusters):
        print(f"\n--- Cluster {i+1} ---")
        c_path = input("Path to this cluster’s .fasta folder? ").strip()
        c_gc_str = input("GC% for this cluster? (e.g. 59.6): ").strip()
        c_name = input("Name of this cluster? (e.g. BD, BE, AK): ").strip()
        try:
            c_gc = float(c_gc_str)
        except ValueError:
            c_gc = None
        cluster_info.append((c_path, c_gc, c_name))

    # -----------------------------
    # 3) Parse the bacterial folder
    # -----------------------------
    print("\nCalculating RSCU for Bacteria...\n")
    bacteria_RSCU = Parse_CDS(bacteria_path)

    # Prepare data array (each entry is a dict of codon -> RSCU)
    # We'll put the bacterial RSCU first
    data = [bacteria_RSCU]
    # Also track names for columns
    col_names = [f"Bacteria (GC={bacteria_gc if bacteria_gc else 'unknown'})"]

    # -----------------------------
    # 4) Parse each cluster
    # -----------------------------
    for (p_path, p_gc, p_name) in cluster_info:
        print(f"\nCalculating RSCU for cluster {p_name} (GC={p_gc})...")
        cluster_RSCU = Parse_CDS(p_path)
        data.append(cluster_RSCU)
        if p_gc:
            col_names.append(f"{p_name} (GC={p_gc})")
        else:
            col_names.append(p_name)

    # -----------------------------
    # 5) Create DataFrame & Heatmap
    # -----------------------------
    print("\nCreating RSCU DataFrame...")
    df = pd.DataFrame.from_dict(data)
    df = df.T  # So codons are rows, and each RSCU dictionary is a column
    df.columns = col_names
    df = df.astype(float)
    print(df)

    # Plot as a heatmap
    plt.figure(figsize=(18, 18))
    plt.imshow(df, interpolation='nearest', cmap='viridis', aspect='auto')
    plt.colorbar()
    plt.yticks(np.arange(len(df.index)), df.index, rotation=0)
    plt.xticks(np.arange(df.shape[1]), df.columns, rotation=90)
    plt.title("RSCU Heatmap")
    plt.savefig("heatmap_RSCU.png")
    plt.close()

    # -----------------------------
    # 6) Compare each cluster vs. bacteria for differences
    # -----------------------------
    difference_data = {}
    # data[0] is the bacterial RSCU
    for i in range(1, len(data)):
        cluster_col_name = col_names[i]
        diff = {}
        for codon in data[0]:
            diff[codon] = abs(data[i][codon] - data[0][codon])
        difference_data[cluster_col_name] = diff

    diff_df = pd.DataFrame.from_dict(difference_data)
    print("\nDifference DataFrame:\n", diff_df)

    # Scatter plot with polynomial fit
    x_array = np.arange(len(diff_df.index))
    colors = list("rgbcmyk")

    plt.figure(figsize=(16, 6))
    plt.rcParams.update({'font.size': 8})

    for i, cluster_col_name in enumerate(diff_df.columns):
        y_array = diff_df[cluster_col_name].values
        color = colors[i % len(colors)]
        # Plot scatter
        plt.scatter(x_array, y_array, label=cluster_col_name, color=color)
        # Fit polynomial (for smoothing)
        coeffs = np.polyfit(x_array, y_array, deg=3)
        poly_fn = np.poly1d(coeffs)
        x_line = np.linspace(x_array.min(), x_array.max(), 100)
        y_line = poly_fn(x_line)
        plt.plot(x_line, y_line, linestyle='-', linewidth=2, color=color)

    plt.xticks(x_array, diff_df.index, rotation=90)
    plt.xlabel("Codon")
    plt.ylabel("RSCU difference vs. Bacteria")
    plt.legend()
    plt.tight_layout()
    plt.savefig("scatterplot_RSCU_difference.png")
    plt.close()

    # If you want to save the frames:
    # df.to_csv('Multi_Cluster_df.csv', index=True)
    diff_df.to_csv('RSCU_Differences_vs_Bacteria.csv', index=True)

    print("\nAll done! Plots and CSVs have been saved.\n")

if __name__ == "__main__":
    main()
