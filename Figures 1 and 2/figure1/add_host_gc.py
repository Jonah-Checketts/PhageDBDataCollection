import pandas as pd

# Load the CSV files
streptomyces_df = pd.read_csv("streptomyces_csv_data.csv")
reproscatter_df = pd.read_csv("phage_host_gc.csv")

# Extract species from strain_name (second and third elements)
streptomyces_df["species"] = streptomyces_df["strain_name"].apply(
    lambda x: " ".join(x.split()[1:3])
)

# Compute average gc_percent for each unique species
avg_gc_by_species = (
    streptomyces_df.groupby("species")["gc_percent"]
    .mean()
    .reset_index()
    .rename(columns={"gc_percent": "avg_gc"})
)

# Add manually retrieved GC content values
manual_gc = pd.DataFrame({
    "species": [
        "Streptomyces azureus",
        "Streptomyces venezuelae",
        "Streptomyces coelicolor"
    ],
    "avg_gc": [0.71, 0.725, 0.72]
})

# Combine automatic and manual GC values
combined_gc = pd.concat([avg_gc_by_species, manual_gc], ignore_index=True)

# Map the average GC content to the REPROSCATTER dataframe
reproscatter_df["gc"] = reproscatter_df["host_species"].map(
    combined_gc.set_index("species")["avg_gc"]
)

# Remove any rows where gc is NaN
reproscatter_df = reproscatter_df.dropna(subset=["gc"])

reproscatter_df["gc"] = reproscatter_df["gc"] * 100

# Save the updated dataframe
reproscatter_df.to_csv("phage_host_gc_full.csv", index=False)
