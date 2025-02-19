from itertools import product
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np

# Load the dataset based on the provided file
phage_host_df = pd.read_csv("heatmap_generation/put in the github/phage_host_data.csv")

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
        gc_diff = None # Handle cases where no match exists
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
plt.figure(figsize=(12, 8))
ax = sns.heatmap(heatmap_data_sorted, cmap='coolwarm_r', annot=False, linewidths=0.5,
                 cbar_kws={'label': 'GC Content Difference'})
plt.title("Phage-Host GC Content Correlation Heatmap (Diagonal Alignment)")
plt.xlabel("Host Species")
plt.ylabel("Phage Name")
plt.xticks(rotation=30, ha="right")
plt.yticks(rotation=0)
plt.tight_layout()

# Highlight known phage-host relationships with black borders in the new diagonal arrangement
for phage, host in known_relationships:
    if phage in heatmap_data_sorted.index and host in heatmap_data_sorted.columns:
        y, x = heatmap_data_sorted.index.get_loc(phage), heatmap_data_sorted.columns.get_loc(host)
        ax.add_patch(plt.Rectangle((x, y), 1, 1, fill=False, edgecolor='black', lw=2))

# Display the updated plot
plt.savefig("phage_host_gc_heatmap.png", dpi=300, bbox_inches='tight')

plt.show()
