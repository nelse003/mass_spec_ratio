import pandas as pd

mounted_df = pd.read_csv("mounted_ratios.csv")
occ_correct_df = pd.read_csv("occ_correct.csv")

occ_df = pd.merge(
    occ_correct_df, mounted_df, right_on="  Mounted Crystal ID ", left_on="crystal"
)
occ_df = occ_df[occ_df.state == "bound"]
occ_df = occ_df[occ_df.resname == "LIG"]
occ_df = occ_df[occ_df["occupancy group"] == 5]
