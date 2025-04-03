import os
import sys
from Bio import SeqIO
from collections import defaultdict
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# Usage: python Calculating_RSCU_Values.py
# Need to edit the script to actually identify clusters and bacteria folder 

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

def clean_cds(filepath): #This is specifically for the bacterial sequences which arent always in multiples of 3.
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
    synonymous_acid = CODON_TABLE[codon]
    for amino_acid in CODON_TABLE:
        if CODON_TABLE[amino_acid] == synonymous_acid:
            count += 1
    return count

def Calculate_RSCU(codon_counts, amino_acid_counts):
    RSCU_data = {}
    #Calculate RSCU    
    #RSCU = (Observed count of codon) / (Total count of codons for the amino acid ÷ Number of synonymous codons)   
    R_val = 0
    for cd in CODON_TABLE:
        R_val = codon_counts[cd] / (amino_acid_counts[CODON_TABLE[cd]] /  Get_Synonymous_Codons(cd))
        #print(f"RSCU val for {cd}: {R_val}")
        RSCU_data[cd] = R_val
    return RSCU_data

def Parse_CDS(cluster_folder):
    #Parsing .fasta files for a given cluster. Each phage CDS is stored in its own .fasta files in the larger folder:
    cluster_folder = os.fsencode(cluster_folder)
    codon_counts = defaultdict(int)
    amino_acid_counts = defaultdict(int)
        
    for file in os.listdir(cluster_folder):
        filename = os.fsdecode(file)
        if not filename.endswith('.fasta'):
            continue
        file_path = os.path.join(os.fsdecode(cluster_folder), filename)
        #clean_cds(file_path)
        for record in SeqIO.parse(file_path, "fasta"):
            try:
                if len(record.seq) % 3 != 0:
                    raise ValueError("Not a multiple of 3 CDS, translation invalid.")
                #print(len(record.seq)/3)
                n = len(record.seq)
                for window in range(0, n, 3):
                    #Count the codons:
                    codon = str(record.seq[window:window+3])
                    if str(codon) in ("TGA", "TAA", "TAG"):
                        continue
                    codon_counts[codon] += 1
                    amino_acid_counts[CODON_TABLE[codon]] += 1

            except ValueError as e:
                print(e)
    dat = Calculate_RSCU(codon_counts, amino_acid_counts)
    return dat

cluster_directories  = [
    "/home/atkizen/Phage_Annotation/CodonUssageAnalyses/cluster_BD_cds", #Infects Streptomyces (exact species) : GC% = 66.4
    "/home/atkizen/Phage_Annotation/CodonUssageAnalyses/cluster_BE_cds", #Infects Streptomyces (mixed) : GC% = 49.6
    "/home/atkizen/Phage_Annotation/CodonUssageAnalyses/cluster_BI_cds", #Infects Streptomyces (mixed) : GC% = 59.6
    "/home/atkizen/Phage_Annotation/CodonUssageAnalyses/cluster_AK_cds", #Infects Arthobacter : GC% = 61.1
    "/home/atkizen/Phage_Annotation/CodonUssageAnalyses/cluster_AS_cds", #Infects Arthobacter : GC% = 66.7
    "/home/atkizen/Phage_Annotation/CodonUssageAnalyses/cluster_CZ_cds", #Infects Gordonia : GC% = 66.4
    "/home/atkizen/Phage_Annotation/CodonUssageAnalyses/cluster_DJ_cds", #Infects Gordonia : GC% = 51.5
    "/home/atkizen/Phage_Annotation/CodonUssageAnalyses/cluster_P__cds", #Infects Mycobacterium : GC% = 67
    "/home/atkizen/Phage_Annotation/CodonUssageAnalyses/cluster_AU_cds", #Infects Arthobacter : GC% = 50.3
    "/home/atkizen/Phage_Annotation/CodonUssageAnalyses/cluster_BU_cds" #Infects Propionibacterium : GC% = 54.1
    
]
bacteria_direct = "/home/atkizen/Phage_Annotation/CodonUssageAnalyses/Bacterial_Sequences" #GC content = ~71%

data = []
bacteria_RSCU = Parse_CDS(bacteria_direct) 
data.append(bacteria_RSCU)

for d in cluster_directories:
    print(f"Calculating RSCU for cluster {d[-6:]}")
    #d is returned as a dictionary
    data.append(Parse_CDS(d))



# Plotting values as heatmap:
df = pd.DataFrame.from_dict(data)
df = df.T
df.columns = ['Streptomyces', 'BD Cluster Phages', 'BE Cluster Phages', 'BI Cluster Phages', 'AK Cluster Phages', 'AS Cluster Phages', 'CZ Cluster Phages', 'DJ Cluster Phages', 'P Cluster Phages', 'AU Cluster Phages', 'BU Cluster phages']
df = df.astype(float)

print(df)

plt.figure(figsize=(18, 18))
plt.imshow(df, interpolation='nearest', cmap='viridis', aspect='auto')
plt.colorbar()
plt.yticks(np.arange(len(df.index)), df.index, rotation=0)
#plt.xticks([0, 1, 2, 3], ['Streptomyces', 'Cluster BD', 'Cluster BE', 'Cluster AK'])  
plt.xticks(np.arange(df.shape[1]), df.columns, rotation=0)
plt.title("RSCU Heatmap")
plt.savefig("heatmatp_RSCU.png")


#Comparing RSCU difference values and plotting on graph:
difference_data = {}
for i in range(1, len(data)):
    diff = {}
    for codon in data[0]:
        diff[codon] = abs(data[i][codon] - data[0][codon])
    difference_data[f"Cluster {cluster_directories[i-1][-6:]}"] = diff
    
diff_df = pd.DataFrame.from_dict(difference_data)
print(diff_df)
colors = list("rgbcmyk")
x_array = np.array(range(len(diff_df.index)))

plt.close('all')
plt.figure(figsize=(16, 6))
plt.rcParams.update({'font.size': 8})

for i, cluster in enumerate(diff_df.columns):
    y_array = diff_df[cluster].values
    color = colors[i % len(colors)]
    # Plot scatter and label only once per cluster.
    plt.scatter(x_array, y_array, label=cluster, color=color)
    # Smoothin a line 
    coeffs = np.polyfit(x_array, y_array, deg=3)
    poly_fn = np.poly1d(coeffs)
    x_line = np.linspace(x_array.min(), x_array.max(), 100)
    y_line = poly_fn(x_line)
    # Plot the fitted line with the same color
    plt.plot(x_line, y_line, linestyle='-', linewidth=2, color=color)

plt.xticks(x_array, diff_df.index, rotation=90)
plt.xlabel("Codon")
plt.ylabel("RSCU difference value")
plt.legend()
plt.savefig("scatterplot_RSCU_difference.png")

#df.to_csv('Multi_Cluster_df.csv', index=True) #If you want just the RSCU values, use this
diff_df.to_csv('RSCU_Difference_From_Bacteria_Multi_Cluster_df.csv', index=True) #If you want a comparison between a bacterial host and the phage clusters, use this