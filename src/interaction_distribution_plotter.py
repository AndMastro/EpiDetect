import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt

df_interactions = pd.DataFrame(columns= ["Method", "SNP rank", "Number of interactions"])
TRAIT = "PP"
plt.grid(visible=True, axis="y", linestyle='-', linewidth=0.5, alpha=0.8)
#####BOOST######
#SBP
# df_interactions.loc[0] = ["BOOST", 1, 226]
# df_interactions.loc[1] = ["BOOST", 2, 180]
# df_interactions.loc[2] = ["BOOST", 3, 157]
# df_interactions.loc[3] = ["BOOST", 4, 130]
# df_interactions.loc[4] = ["BOOST", 5, 89]

#DBP
# df_interactions.loc[0] = ["BOOST", 1, 323]
# df_interactions.loc[1] = ["BOOST", 2, 247]
# df_interactions.loc[2] = ["BOOST", 3, 183]
# df_interactions.loc[3] = ["BOOST", 4, 120]
# df_interactions.loc[4] = ["BOOST", 5, 98]

#PP
df_interactions.loc[0] = ["BOOST", 1, 212]
df_interactions.loc[1] = ["BOOST", 2, 160]
df_interactions.loc[2] = ["BOOST", 3, 149]
df_interactions.loc[3] = ["BOOST", 4, 139]
df_interactions.loc[4] = ["BOOST", 5, 92]

#####MDR######
#SBP
# df_interactions.loc[5] = ["MDR", 1, 263]
# df_interactions.loc[6] = ["MDR", 2, 263]
# df_interactions.loc[7] = ["MDR", 3, 237]
# df_interactions.loc[8] = ["MDR", 4, 163]
# df_interactions.loc[9] = ["MDR", 5, 46]

#DBP
# df_interactions.loc[5] = ["MDR", 1, 252]
# df_interactions.loc[6] = ["MDR", 2, 228]
# df_interactions.loc[7] = ["MDR", 3, 181]
# df_interactions.loc[8] = ["MDR", 4, 139]
# df_interactions.loc[9] = ["MDR", 5, 72]

#PP
df_interactions.loc[5] = ["MDR", 1, 282]
df_interactions.loc[6] = ["MDR", 2, 282]
df_interactions.loc[7] = ["MDR", 3, 282]
df_interactions.loc[8] = ["MDR", 4, 151]
df_interactions.loc[9] = ["MDR", 5, 10]

#####NID######
#SBP
# df_interactions.loc[10] = ["NID", 1, 163]
# df_interactions.loc[11] = ["NID", 2, 115]
# df_interactions.loc[12] = ["NID", 3, 114]
# df_interactions.loc[13] = ["NID", 4, 92]
# df_interactions.loc[14] = ["NID", 5, 91]

#DBP
# df_interactions.loc[10] = ["NID", 1, 161]
# df_interactions.loc[11] = ["NID", 2, 122]
# df_interactions.loc[12] = ["NID", 3, 102]
# df_interactions.loc[13] = ["NID", 4, 87]
# df_interactions.loc[14] = ["NID", 5, 77]

#PP
df_interactions.loc[10] = ["NID", 1, 178]
df_interactions.loc[11] = ["NID", 2, 98]
df_interactions.loc[12] = ["NID", 3, 63]
df_interactions.loc[13] = ["NID", 4, 63]
df_interactions.loc[14] = ["NID", 5, 56]

#####EpiCID######
#SBP
# df_interactions.loc[15] = ["EpiCID", 1, 150]
# df_interactions.loc[16] = ["EpiCID", 2, 118]
# df_interactions.loc[17] = ["EpiCID", 3, 68]
# df_interactions.loc[18] = ["EpiCID", 4, 66]
# df_interactions.loc[19] = ["EpiCID", 5, 64]

#DBP
# df_interactions.loc[15] = ["EpiCID", 1, 121]
# df_interactions.loc[16] = ["EpiCID", 2, 79]
# df_interactions.loc[17] = ["EpiCID", 3, 76]
# df_interactions.loc[18] = ["EpiCID", 4, 70]
# df_interactions.loc[19] = ["EpiCID", 5, 69]

#PP
df_interactions.loc[15] = ["EpiCID", 1, 130]
df_interactions.loc[16] = ["EpiCID", 2, 103]
df_interactions.loc[17] = ["EpiCID", 3, 85]
df_interactions.loc[18] = ["EpiCID", 4, 60]
df_interactions.loc[19] = ["EpiCID", 5, 59]


sns.barplot(x = df_interactions["SNP rank"], y=df_interactions["Number of interactions"], hue = df_interactions["Method"])
plt.savefig("../data/results/epistaticInteraction/comparison_plots/" + TRAIT + "/interactions_distribution_thesis_28_10_2023.pdf", dpi=300)
plt.savefig("../data/results/epistaticInteraction/comparison_plots/" + TRAIT + "/interactions_distribution_thesis_28_10_2023.eps", dpi=300)
plt.savefig("../data/results/epistaticInteraction/comparison_plots/" + TRAIT + "/interactions_distribution_thesis_28_10_2023.png", dpi=300)

# sns.kdeplot(x = df_interactions["SNP RANK"], y=df_interactions["# interactions"], hue = df_interactions["Method"])
# plt.savefig("interactions_distribution_kde.png")


print(df_interactions)