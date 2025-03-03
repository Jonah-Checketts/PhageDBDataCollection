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
    x='cluster',
    y='GC_percent',
    palette='Set1',       # Choose a color palette
    dodge=True,           # Separate the points by hue in each category
    jitter=True,
    size=5
)

# Label the axes
plt.xlabel("Cluster")
plt.ylabel("GC Percent")

# Optionally, give the plot a title
plt.title("GC Percent by Cluster")

# Display the plot
plt.show()
