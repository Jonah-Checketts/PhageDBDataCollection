#!/usr/bin/env python3

import re
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.colors import Normalize

def main():
    # -------------------------------------------------------------------
    # 1) Prompt the user
    # -------------------------------------------------------------------
    host_gc_str = input("Enter the host (bacteria) GC% (e.g. '71'): ").strip()
    try:
        host_gc = float(host_gc_str)
    except ValueError:
        print("Invalid GC% entered. Defaulting to 71.0")
        host_gc = 71.0

    infect_str = input(
        "Enter a comma-separated list of clusters that infect this host\n"
        "(e.g. 'AK, AS, BD'): "
    ).strip()
    infecting_clusters = {s.strip() for s in infect_str.split(",") if s.strip()}

    # -------------------------------------------------------------------
    # 2) Read the CSV
    # -------------------------------------------------------------------
    csv_file = "RSCU_Differences_vs_Bacteria.csv"
    print(f"Reading difference data from {csv_file}...")
    data = pd.read_csv(csv_file, index_col=0)

    # -------------------------------------------------------------------
    # 3) Parse each column name for cluster name + GC
    #    Example: "AK (GC=61.1)"
    # -------------------------------------------------------------------
    pattern = re.compile(r"^(.+?)\s*\(GC\s*=\s*([\d\.]+)\)$")

    cluster_names = {}
    cluster_gcs = {}

    for col in data.columns:
        match = pattern.match(col)
        if match:
            c_name = match.group(1).strip()    
            c_gc   = float(match.group(2))     
        else:
            c_name = col
            c_gc   = 0.0
        cluster_names[col] = c_name
        cluster_gcs[col]   = c_gc

    # -------------------------------------------------------------------
    # 4) Color each cluster by difference from host GC
    #    Now from 2.5 to 25
    # -------------------------------------------------------------------
    norm = Normalize(vmin=2.5, vmax=25)
    cmap = plt.get_cmap("viridis_r")

    def cluster_color(col_name):
        diff_gc = abs(cluster_gcs[col_name] - host_gc)
        return cmap(norm(diff_gc))

    # -------------------------------------------------------------------
    # 5) Plot
    # -------------------------------------------------------------------
    x_array = np.arange(len(data.index))
    plt.figure(figsize=(10, 6))
    plt.rcParams.update({'font.size': 9})

    for col in data.columns:
        y_array = data[col].values
        color   = cluster_color(col)

        # Solid if in infecting_clusters, else dashed
        c_name = cluster_names[col]
        linestyle = '-' if c_name in infecting_clusters else ':'

        plt.scatter(x_array, y_array, color=color, s=10)
        coeffs = np.polyfit(x_array, y_array, deg=3)
        poly_fn = np.poly1d(coeffs)
        x_line = np.linspace(x_array.min(), x_array.max(), 100)
        y_line = poly_fn(x_line)
        plt.plot(x_line, y_line, linestyle=linestyle, linewidth=2, color=color, label=col)

    plt.xticks(x_array, data.index, rotation=90, fontsize=6)
    plt.xlabel("Codon")
    plt.ylabel("RSCU Difference vs. Reference Bacteria")
    plt.title("RSCU Differences Across Clusters")

    # -------------------------------------------------------------------
    # 6) Legends
    # -------------------------------------------------------------------
    cluster_legend = plt.legend(title="Clusters", loc="upper left", fontsize=7)
    plt.gca().add_artist(cluster_legend)

    solid_handle = Line2D([0], [0], color='gray', lw=2, linestyle='-',
                          label='Solid: Infects Host')
    dashed_handle = Line2D([0], [0], color='gray', lw=2, linestyle=':',
                           label="Dashed: Doesn't Infect Host")
    plt.legend(handles=[solid_handle, dashed_handle],
               title="Line Style Key", loc="upper right")

    # -------------------------------------------------------------------
    # 7) Colorbar
    # -------------------------------------------------------------------
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    cbar = plt.colorbar(sm, ax=plt.gca(), pad=0.02)
    cbar.set_label("GC% Difference vs. Host")

    plt.tight_layout()
    plt.savefig("scatterplot_RSCU_difference.png")
    plt.close()

    print("Plot saved to scatterplot_RSCU_difference.png")

if __name__ == "__main__":
    main()
