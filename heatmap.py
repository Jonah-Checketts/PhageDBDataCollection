import matplotlib.pyplot as plt
import numpy as np
import matplotlib
import matplotlib as mpl
import math
from matplotlib.colors import LinearSegmentedColormap

num_bac = 20
num_phage = 20
num_phage_cat = 1
num_bac_cat = 1

phage_info = {}
bacteria_info = {}

first_line = True
with open("cluster_BD_data.csv", "r") as in_f:
    i = 0
    for line in in_f:
        if first_line:
            first_line = False
            continue
        if i > num_phage:
            break
        i += 1
        cols = line.strip().split(",")
        phage_info[cols[0]] = cols[2]

'''first_line = True
with open("cluster_AW_data.csv", "r") as in_f:
    i = 0
    for line in in_f:
        if first_line:
            first_line = False
            continue
        if i > num_phage:
            break
        i += 1
        cols = line.strip().split(",")
        phage_info[cols[0]] = cols[2]'''

first_line = True
with open("streptomyces_csv_data.csv") as in_f:
    i = 0
    for line in in_f:
        if first_line:
            first_line = False
            continue
        if i > num_bac:
            break
        i += 1
        cols = line.strip().split(",")
        bacteria_info[cols[0]] = round((float(cols[1]) * 100), 2)

bacteria_names = list(bacteria_info.keys())
phage_names = list(phage_info.keys())

heatmap_data = []

for bacteria_name, bacteria_gc_per in bacteria_info.items():
    curr_row = []
    for phage_name, phage_gc_per in phage_info.items():
        curr_row.append(abs(float(bacteria_gc_per) - float(phage_gc_per)))
    heatmap_data.append(curr_row)

heatmap_data = np.array(heatmap_data)
fig, ax = plt.subplots()
im = ax.imshow(heatmap_data, cmap="RdBu_r")

ax.set_yticks(range(num_bac * num_bac_cat), labels=bacteria_names[:(num_bac * num_bac_cat)])
ax.set_xticks(range(num_phage * num_phage_cat), labels=phage_names[:(num_phage * num_phage_cat)], rotation=90, ha="right", rotation_mode="anchor")

cbar = ax.figure.colorbar(im, ax=ax)
cbar.ax.set_ylabel("Difference in phage-host GC percent", rotation=-90, va="bottom")

ax.set_title("Phage-Host GC Percent Coorelation")
fig.tight_layout()
plt.show()