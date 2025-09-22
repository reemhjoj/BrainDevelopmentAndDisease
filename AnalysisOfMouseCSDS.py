import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import mannwhitneyu

# Load dataset (downloaded from Gene expression omnibus GSE146845), data file was manually filled based on the data
# after creating a csv file of the data it was uploaded for plotting and analysis

df_mouse = pd.read_csv(r"C:\Users\...\TMEM132D_mouse_GSE146845.csv")



# filtering the groups to include control and chronic stress model only
df_filtered = df_mouse[df_mouse["Group"].isin(["CTRL", "stress/vehicle "])].copy()

# Rename chronic stress model group for clearer plotting
df_filtered["Group"] = df_filtered["Group"].replace({"stress/vehicle ": "CSDS"})

# statistics with Mann-Whitney U test
ctrl_vals = df_filtered[df_filtered["Group"] == "CTRL"]["Expression"]
csds_vals = df_filtered[df_filtered["Group"] == "CSDS"]["Expression"]
stat, pval = mannwhitneyu(ctrl_vals, csds_vals, alternative="two-sided")
print(stat, pval)
# Determine significance stars for plotting
if pval < 0.001:
    star = "***"
elif pval < 0.01:
    star = "**"
elif pval < 0.05:
    star = "*"
else:
    star = "ns"

# Plotting boxplot for CTRL and CSDS groups
plt.figure(figsize=(8,6))
ax = sns.boxplot(data=df_filtered, x="Group", y="Expression",
                 order=["CTRL", "CSDS"], palette=["white", "lightgray"],
                 showfliers=False, linewidth=1.5)

# plotting datapoints of each subject in each group
sns.stripplot(data=df_filtered, x="Group", y="Expression",
              order=["CTRL", "CSDS"], color="black", size=5, jitter=True)


# adjusting the location of the stars and the bar
group_max = df_filtered.groupby("Group")["Expression"].max()
y_max = group_max.max()
y_min = df_filtered["Expression"].min()
bar_height = y_max + 0.05
text_height = y_max + 0.10
x1, x2 = 0, 1
ax.plot([x1, x1, x2, x2], [bar_height, bar_height, bar_height, bar_height], color="black", linewidth=1.2)
ax.text((x1+x2)/2, text_height, star, ha="center", va="bottom", fontsize=14, weight="bold")

# adjusting the size of ticks at x-axis and y-axis
ax.tick_params(axis='both', which='major', labelsize=18)
# Adding extra space on x and y-axis for clearer figure
ax.set_ylim(y_min - 0.1, y_max + 0.3)
# adding titles
plt.ylabel("Tmem132d expression (normalized CPM)",  fontsize=20)
plt.xlabel("")
plt.title("PFC",  fontsize=16, weight="bold")
#saving the plot
plt.savefig(r"C:\Users\...", dpi=300, bbox_inches="tight")
plt.tight_layout()
plt.show()