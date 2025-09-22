import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import mannwhitneyu

# Loading human data file (GSE80655) which was downloaded from Gene expression ominbus
# the datafile analyzed here was manually created based on the downloaded file and GEO website.
df_human = pd.read_csv(r"C:\Users\...")
# calculting log2 of FPKM for plotting
df_human["log2_FPKM"] = np.log2(df_human["TMEM132D_FPKM"])
# our Regions of interest are the anterior cingulate gyrus (AnCg) and the dorsal lateral prefrontal cortex (DLPFC)
# these two ROIs are filtered from the database
df_human = df_human[df_human["Brain region"].isin(["DLPFC", "AnCg"])]
# Lableing adjusment for the plotting
label_map = {"CTRL": "CTRL", "MDD": "MDD", "Schizophrenia": "SCZ", "Bipolar": "BPD"}
df_human["Group_short"] = df_human["Group"].map(label_map)
# creating a unified group order for the plotting
group_order = ["CTRL", "MDD", "SCZ", "BPD"]


# This function takes the data file and a region of interest and plots the data as a boxplot
def plot_region(ax, df, region_name):
   # filtering the region of interst
    sub = df[df["Brain region"] == region_name]

    # plotting a boxplot to visualize the different psychiatric disorders
    sns.boxplot(
        data=sub, x="Group_short", y="log2_FPKM",
        order=group_order, palette=["white"] * 4,
        showcaps=True, fliersize=0,
        linewidth=1.8, ax=ax
    )

    # Keeping the plot black and white (no colors)
    for i, patch in enumerate(ax.artists):
        patch.set_facecolor("white")
        patch.set_edgecolor("black")
        patch.set_linewidth(1.8)

    # plotting the log2 FPKM datapoints values on the box plot for different indivduals
    sns.stripplot(data=sub, x="Group_short", y="log2_FPKM",order=group_order, dodge=True,color="black", size=4, jitter=True, ax=ax)

    # Statistical testing : Mann-Whitney tests vs CTRL
    ctrl_vals = sub[sub["Group_short"] == "CTRL"]["log2_FPKM"]
    for disorder in ["MDD", "SCZ", "BPD"]:
        dis_vals = sub[sub["Group_short"] == disorder]["log2_FPKM"]
        if len(dis_vals) > 0:
            stat, pval = mannwhitneyu(ctrl_vals, dis_vals, alternative="two-sided")
            print(disorder)
            print(stat, pval)
            if pval < 0.05:
                if pval < 0.001:
                    star = "***"
                elif pval < 0.01:
                    star = "**"
                else:
                    star = "*"
                ymax = max(sub["log2_FPKM"]) + 0.3
                xpos = group_order.index(disorder)
                ax.text(xpos, ymax, star, ha="center", va="bottom", fontsize=14, weight="bold")

    ax.set_title(region_name, fontsize=16, weight="bold")
    ax.set_xlabel("")
    ax.tick_params(axis='both', labelsize=18)


# plotting the two regions plots in one figure
fig, axes = plt.subplots(1, 2, figsize=(8, 6), sharey=True)

plot_region(axes[0], df_human, "AnCg")
plot_region(axes[1], df_human, "DLPFC")

# only one y-axis title
axes[0].set_ylabel("TMEM132D expression (log2 FPKM)", fontsize=20)
axes[1].set_ylabel("")


plt.tight_layout(rect=[0, 0, 1, 0.95])
# saving the plot
plt.savefig(r"C:\Users\...", dpi=300, bbox_inches="tight")

plt.show()
