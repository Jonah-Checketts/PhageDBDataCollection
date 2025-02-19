from itertools import product
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.patches as mpatches

# Load the dataset based on the provided file
phage_host_df = pd.read_csv("/Users/blakemcgee/Desktop/heatmap/phage_host_data.csv")

# Inspect the dataset to ensure column names are correct
print("Dataset Columns:", phage_host_df.columns)

# Extract unique phages and hosts
phages = phage_host_df['Phage Name'].unique()
hosts = phage_host_df['Species'].unique()

# Create a full phage-host GC content difference matrix
full_gc_data = []
known_relationships = set(zip(phage_host_df['Phage Name'], phage_host_df['Species']))

for phage, host in product(phages, hosts):
    phage_gc = phage_host_df.loc[phage_host_df['Phage Name'] == phage, 'Phage GC'].values
    host_gc = phage_host_df.loc[phage_host_df['Species'] == host, 'Host GC'].values
    if len(phage_gc) > 0 and len(host_gc) > 0:
        gc_diff = abs(phage_gc[0] - host_gc[0])
    else:
        gc_diff = None  # Handle cases where no match exists
    full_gc_data.append([phage, host, gc_diff])

# Convert to DataFrame
full_gc_df = pd.DataFrame(full_gc_data, columns=['Phage Name', 'Species', 'GC_Content_Difference'])

# Pivot table for heatmap
heatmap_data = full_gc_df.pivot(index='Phage Name', columns='Species', values='GC_Content_Difference')

# Reorder phages and hosts to align black boxes along a diagonal
phage_order = full_gc_df.groupby("Phage Name")["GC_Content_Difference"].mean().sort_values().index
host_order = full_gc_df.groupby("Species")["GC_Content_Difference"].mean().sort_values().index

# Reorder the heatmap data
heatmap_data_sorted = heatmap_data.loc[phage_order, host_order]

# Plot the updated heatmap
plt.figure(figsize=(12, 12))  # Increase width and height

ax = sns.heatmap(heatmap_data_sorted, cmap='coolwarm_r', annot=False, linewidths=0.5,
                 cbar_kws={'label': 'GC Content Difference'})

plt.title("Phage-Host GC Content Correlation Heatmap (Diagonal Alignment)")
plt.xlabel("Host Species")
plt.ylabel("Phage Name")
plt.xticks(rotation=30, ha="right")
plt.yticks(rotation=0)

# Adjust layout to make room for the text
plt.subplots_adjust(bottom=0.25)  # Allocate space at the bottom

# Highlight known phage-host relationships with black borders
for phage, host in known_relationships:
    if phage in heatmap_data_sorted.index and host in heatmap_data_sorted.columns:
        y, x = heatmap_data_sorted.index.get_loc(phage), heatmap_data_sorted.columns.get_loc(host)
        ax.add_patch(plt.Rectangle((x, y), 1, 1, fill=False, edgecolor='black', lw=2))

# Add a longer figure description below
description = ("Figure 1: This heatmap visualizes the GC content differences between phages and their potential host species. "
               "Darker red colors indicate a lesser GC content difference, while darker blue shades represent larger differences. "
               "Black-bordered cells highlight known phage-host relationships. "
               "This visualization shows that within the BD cluster of phages, the GC content difference has no correlation with host specificity. "
               "The claim that GC content similarity is a predictor of a phage-host relationship is not supported by this data of BD cluster phages. "
               "https://github.com/Jonah-Checketts/PhageDBDataCollection/tree/a27d9d0fcc361a8571811fcf655ee3b3128dcda0/heatmap_generation/put%20in%20the%20github")

plt.figtext(0.5, 0.005, description, wrap=True, horizontalalignment='center', fontsize=10, bbox=dict(facecolor='white', alpha=0.5))

# Save with additional padding
plt.savefig("phage_host_gc_heatmap.png", dpi=300, bbox_inches='tight', pad_inches=1)

plt.show()