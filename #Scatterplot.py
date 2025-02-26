#Scatterplot

import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd

df_updated_bdhost = pd.read_csv("Updated_BDHost_Data__GC_Multiplied_by_100_.csv")

# Drop rows with missing GC values
df_cleaned = df_updated_bdhost.dropna(subset=['GC_percent', 'gc'])

# Convert columns to numeric if they are not already
df_cleaned['GC_percent'] = pd.to_numeric(df_cleaned['GC_percent'], errors='coerce')
df_cleaned['gc'] = pd.to_numeric(df_cleaned['gc'], errors='coerce')

# Create scatter plot
plt.figure(figsize=(8, 6))
sns.scatterplot(x=df_cleaned['GC_percent'], y=df_cleaned['gc'], alpha=0.7, label="Phage-Host GC%")

# Fit and plot regression line
sns.regplot(x=df_cleaned['GC_percent'], y=df_cleaned['gc'], scatter=False, color='red', label="Trend Line")

# Labels and title
plt.xlabel("Phage GC%")
plt.ylabel("Host GC%")
plt.title("Scatter Plot of Phage GC% vs. Host GC%")
plt.legend()
plt.grid(True)

# Show the plot
plt.show()
