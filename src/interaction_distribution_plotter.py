import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt

df_interactions = pd.DataFrame(columns= ["Method", "SNP RANK", "# interactions"])

#####BOOST######
df_interactions.loc[0] = ["BOOST", 1, 226]
df_interactions.loc[1] = ["BOOST", 2, 180]
df_interactions.loc[2] = ["BOOST", 3, 157]
df_interactions.loc[3] = ["BOOST", 4, 130]
df_interactions.loc[4] = ["BOOST", 5, 89]

#####MDR######
df_interactions.loc[5] = ["MDR", 1, 263]
df_interactions.loc[6] = ["MDR", 2, 263]
df_interactions.loc[7] = ["MDR", 3, 237]
df_interactions.loc[8] = ["MDR", 4, 163]
df_interactions.loc[9] = ["MDR", 5, 46]

#####NID######
df_interactions.loc[10] = ["NID", 1, 193]
df_interactions.loc[11] = ["NID", 2, 101]
df_interactions.loc[12] = ["NID", 3, 98]
df_interactions.loc[13] = ["NID", 4, 80]
df_interactions.loc[14] = ["NID", 5, 72]

#####MDR######
df_interactions.loc[15] = ["EpiCID", 1, 142]
df_interactions.loc[16] = ["EpiCID", 2, 125]
df_interactions.loc[17] = ["EpiCID", 3, 92]
df_interactions.loc[18] = ["EpiCID", 4, 75]
df_interactions.loc[19] = ["EpiCID", 5, 72]



sns.barplot(x = df_interactions["SNP RANK"], y=df_interactions["# interactions"], hue = df_interactions["Method"])
plt.savefig("interactions_distribution_300dpi.png", dpi=300)

# sns.kdeplot(x = df_interactions["SNP RANK"], y=df_interactions["# interactions"], hue = df_interactions["Method"])
# plt.savefig("interactions_distribution_kde.png")


print(df_interactions)