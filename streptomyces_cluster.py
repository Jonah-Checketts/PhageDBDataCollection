import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# Read in your data
df = pd.read_csv("streptomyces_phages.csv")

# Set a clean style
sns.set_style("whitegrid")

# Create a strip plot (jittered points) with GC_percent on the y-axis
# and cluster categories on the x-axis
sns.stripplot(
    data=df,
    y='GC_percent',      # Choose a color palette
    dodge=True,           # Separate the points by hue in each category
    jitter=True,
    size=5
)



# Label the axes
plt.xlabel("Phage")
plt.ylabel("GC Percent")

# Optionally, give the plot a title
plt.title("GC Percent of Streptomyces Phages")

# Display the plot
plt.show()
