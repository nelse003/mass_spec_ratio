import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np

def closest_peak(peak_val):
    peaks = [25124, 24218, 23525, 22786]

    return min(peaks, key=lambda x:abs(x-peak_val))


post_crystal_df = pd.read_csv("post_crystal_data_plate_well.csv")

post_crystal_df.drop('calibrated_ratio', inplace=True, axis=1)
post_crystal_df.drop('Unnamed: 0', inplace=True, axis=1)
post_crystal_df.drop('Unnamed: 0.1', inplace=True, axis=1)
post_crystal_df.drop_duplicates(inplace=True)

#post_crystal_df = post_crystal_df[~post_crystal_df["key"].str.contains("NU0000308a")]
#post_crystal_df = post_crystal_df[~post_crystal_df["key"].str.contains("mg_ml")]
#post_crystal_df = post_crystal_df[~post_crystal_df["key"].str.contains("No_cov")]
#post_crystal_df = post_crystal_df.rename(index=str,columns={'Unnamed: 0':"peak"})

print(post_crystal_df['peak_loc'].unique())

post_crystal_df['peak'] = post_crystal_df['peak_loc'].apply(closest_peak)

post_crystal_df["ratio_diff"] = abs(post_crystal_df["height_ratio"] - post_crystal_df["intended_ratio"])

std_height_ratios = []
std_area_ratios = []
for key in post_crystal_df['key']:
    key_df =post_crystal_df[post_crystal_df['key']==key]
    std_height_ratios.append(key_df['height_ratio'].std())
    std_area_ratios.append(key_df['area_ratio'].std())

print(np.mean(std_height_ratios), np.mean(std_area_ratios))

ax = plt.subplot(111)

ax.spines["right"].set_visible(False)
ax.spines["top"].set_visible(False)

plt.hist(std_height_ratios, bins=20)
plt.xlabel("Standard deviation between height ratios for each dataset",fontsize=14)
plt.ylabel("Frequency",fontsize=14)
plt.savefig("std_height_ratios",dpi=300)
plt.close()

ax = plt.subplot(111)

ax.spines["right"].set_visible(False)
ax.spines["top"].set_visible(False)

plt.hist(std_area_ratios, bins=20)
plt.xlabel("Standard deviation between area ratios for each dataset",fontsize=14)
plt.ylabel("Frequency",fontsize=14)
plt.savefig("std_area_ratios",dpi=300)
plt.close()


# fig = plt.figure(figsize=(8, 4.5))
#
# ax = sns.violinplot(x="peak", y="ratio_diff", data=post_crystal_df)
#
# post_crystal_df.to_csv("post_crystal_ratio_diff.csv")
# plt.xlabel("Unlabelled peak location (Daltons)")
# plt.ylabel("Difference between Intended ratio and ")
#
# plt.show()