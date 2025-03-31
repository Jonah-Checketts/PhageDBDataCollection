import os
import pandas as pd
import numpy as np
from matplotlib.colors import LinearSegmentedColormap, Normalize
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

data = pd.read_csv('Difference_From_Streptomyces_Multi_Cluster_df.csv', index_col=0)
data.columns = [col.replace("_cds", "") for col in data.columns]
data.rename(columns={'Cluster P_': 'Cluster P'}, inplace=True)
#print(data)

'''
#Scatterplot grouping "infects streptomyces" vs "doesnt infect streptomyces"
x_array = np.array(range(len(data.index)))
plt.figure(figsize=(16, 6))
plt.rcParams.update({'font.size': 8})

for i, cluster in enumerate(data.columns):
    #print(cluster)
    y_array = data[cluster].values
    if cluster in ("Cluster BD_cds", "Cluster BE_cds", "Cluster BI_cds"):
        color = "r"
    else:
        color = "b"
    #color = colors['red']
    # Plot scatter and label only once per cluster.
    #plt.scatter(x_array, y_array, label=cluster, color=color)
    # Smoothin a line 
    coeffs = np.polyfit(x_array, y_array, deg=3)
    poly_fn = np.poly1d(coeffs)
    x_line = np.linspace(x_array.min(), x_array.max(), 100)
    y_line = poly_fn(x_line)
    # Plot the fitted line with the same color. Comment it out to remove smoothed lines. To ONLY have smooth lines, comment out scatter above, and add "label = cluster" below
    plt.plot(x_line, y_line, linestyle='-', linewidth=2, color=color, label = cluster)

plt.xticks(ticks=x_array, labels=data.index, rotation=90, ha='center')
plt.xlabel("Codon")
plt.ylabel("RSCU difference value")
plt.legend()
plt.savefig("scatterplot_RSCU_difference.png")
'''

#Scatterplot heatmap for GC content
#BD GC% = 66.4
#BE GC% = 49.6
#BI GC% = 59.6
#AK GC% = 61.1
#AS GC% = 66.7
#CZ GC% = 66.4
#DJ GC% = 51.5
#P GC% = 67
GC_dict = {"Cluster BD" : abs(66.4 - 71), "Cluster BE" : abs(49.6 - 71), "Cluster BI" : abs(59.6 - 71), "Cluster AK" : abs(61.1 - 71), 
    "Cluster AS" : abs(66.7 - 71), "Cluster CZ" : abs(66.4 - 71), "Cluster DJ" : abs(51.5 - 71), "Cluster P" : abs(67 - 71), "Cluster AU" : abs(50.3 - 71), "Cluster BU" : abs(54.1 - 71)}
norm = Normalize(vmin=2.5, vmax=25)
cmap = plt.get_cmap("viridis_r")
colors = [cmap(norm(v)) for v in GC_dict.values()]
#GC Content Testing
#plt.figure(figsize=(10, 6))
#plt.scatter(range(len(GC_dict)), list(GC_dict.values()), color=colors, s=100)
#plt.xticks(range(len(GC_dict)), list(GC_dict.keys()), rotation=45, ha='right')
#plt.xlabel("Cluster")
#plt.ylabel("GC Content")
#plt.title("GC Content by Cluster")
#lt.savefig("scatterplot_RSCU_difference.png")

new_order = [
    "Cluster BE",  # col1
    "Cluster BU",  # col2
    "Cluster AU",  # col3
    "Cluster DJ",  # col4
    "Cluster CZ",  # col5
    "Cluster BI",  # col6
    "Cluster AK",  # col7
    "Cluster AS",  # col8
    "Cluster P",  # col9
    "Cluster BD"   # col10
]
# Reorder the DataFrame columns
data = data[new_order]

x_array = np.array(range(len(data.index)))

plt.figure(figsize=(8.2, 5.7))
plt.rcParams.update({'font.size': 10})

dashed_clusters = {"Cluster BU", "Cluster AU", "Cluster DJ", "Cluster CZ", "Cluster AK", "Cluster AS", "Cluster P"}

for i, cluster in enumerate(data.columns):
    #print(cluster)
    y_array = data[cluster].values
    gc_value = GC_dict.get(cluster, None)
    if gc_value is not None:
         color = cmap(norm(gc_value))

    #color = colors['red']
    # Plot scatter and label only once per cluster.
    plt.scatter(x_array, y_array, color=color, s = 12)
    # Smoothin a line 
    coeffs = np.polyfit(x_array, y_array, deg=3)
    poly_fn = np.poly1d(coeffs)
    x_line = np.linspace(x_array.min(), x_array.max(), 100)
    y_line = poly_fn(x_line)

    linestyle = ':' if cluster in dashed_clusters else '-'
    # Plot the fitted line with the same color. Comment it out to remove smoothed lines. To ONLY have smooth lines, comment out scatter above, and add "label = cluster" below
    plt.plot(x_line, y_line, linestyle=linestyle, linewidth=2, color=color, label = cluster)

plt.xticks(ticks=x_array, labels=data.index, rotation=90, ha='center', fontsize = 6)
plt.xlabel("Codon")
plt.ylabel("RSCU difference from Streptomyces")
plt.legend()

cluster_legend = plt.legend(title="Clusters", loc="upper left", fontsize = 7)
plt.gca().add_artist(cluster_legend)  # keep the cluster legend on the plot

solid_handle = Line2D([0], [0], color='gray', lw=2, linestyle='-', label='Solid Lines: Infects Streptomyces')
dashed_handle = Line2D([0], [0], color='gray', lw=2, linestyle=':', label='Dashed lines: Does NOT infet Streptomyces')
plt.legend(handles=[solid_handle, dashed_handle], title="Line Style Key", loc="upper right")

plt.tight_layout(pad=1.0)  # padding adjustment

ax = plt.gca()  # get current axes
sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
sm.set_array([])  # required for the colorbar to work
cbar = plt.colorbar(sm, ax=ax, pad=-.002)
cbar.set_label("GC% Difference From Streptomyces")

plt.savefig("scatterplot_RSCU_difference.png")
