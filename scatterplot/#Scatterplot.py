#Scatterplot

import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import statsmodels.api as sm

df_updated_bdhost = pd.read_csv("scatterplot/Updated_BDHost_Data__GC_Multiplied_by_100_.csv")

# Drop rows with missing GC values
df_cleaned = df_updated_bdhost.dropna(subset=['GC_percent', 'gc'])

# Convert columns to numeric if they are not already
df_cleaned['GC_percent'] = pd.to_numeric(df_cleaned['GC_percent'], errors='coerce')
df_cleaned['gc'] = pd.to_numeric(df_cleaned['gc'], errors='coerce')

# Create scatter plot
plt.figure(figsize=(10, 8))
sns.scatterplot(x=df_cleaned['GC_percent'], y=df_cleaned['gc'], alpha=0.7, label="Phage-Host GC%")

# Fit and plot regression line
sns.regplot(x=df_cleaned['GC_percent'], y=df_cleaned['gc'], scatter=False, color='red', label="Trend Line")

# Calculate R-squared value
X = sm.add_constant(df_cleaned['GC_percent'])
model = sm.OLS(df_cleaned['gc'], X).fit()
r_squared = model.rsquared

# Labels and title
plt.xlabel("Phage GC%")
plt.ylabel("Host GC%")
plt.title(f"BD Cluster Phage-Host GC%")

# Add a detailed legend
legend_text = (
    "This scatter plot shows the relationship between the GC% of phages within the BD cluster and their hosts. "
    "Each point represents a phage-host pair, with the x-axis representing the GC% of the phage "
    "and the y-axis representing the GC% of the host. The red trend line indicates the lack "
    "of correlation between phage GC% and host GC% within the members of the BD cluster."
)
plt.figtext(0.5, 0.01, legend_text, wrap=True, horizontalalignment='center', fontsize=10,
            bbox=dict(boxstyle='round,pad=0.5', edgecolor='black', facecolor='white'))

# Add R-squared value next to the legend box
plt.text(0.95, 0.925, f'R² = {r_squared:.2f}', transform=plt.gca().transAxes, fontsize=10,
         verticalalignment='top', horizontalalignment='right', bbox=dict(boxstyle='round,pad=0.5', edgecolor='black', facecolor='white'))

plt.grid(True)

# Show the plot
plt.show()