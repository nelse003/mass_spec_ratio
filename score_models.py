import os
import pandas as pd

root = "/dls/science/groups/i04-1/elliot-dev/Work/NUDT7A_mass_spec_refinements/"

# Define the file name and the folders to loop over
folders = {
    "copy_atoms/phenix/2019-06-01": ("refine_001.pdb", "refine_001.mtz", "phenix"),
    "copy_atoms/refmac/2019-07-08": ("refine.pdb", "refine.mtz", "refmac"),
    "copy_atoms/buster/2019-05-29": ("refine.pdb", "refine.mtz", "buster"),
    "copy_atoms_190525_buster": ("refine.pdb", "refine.mtz", "buster_superposed"),
    "copy_atoms/refmac_superposed/190706": (
        "refine.pdb",
        "refine.mtz",
        "refmac_superposed",
    ),
    "copy_atoms_190525_phenix": ("refine.pdb", "refine.mtz", "phenix_superposed"),
    "copy_atoms/phenix_b_fix/2019-10-05": ("refine.pdb", "refine.mtz","phenix_b_fix"),
    "copy_atoms/phenix_b_fix_non_superposed/2019-10-04":("refine_001.pdb",
                                                         "refine_001.mtz",
                                                         "phenix_b_fix_non_superposed")

}


scores_df_list = []
for folder in folders:
    for crystal_folder in os.listdir(os.path.join(root, folder)):
        if os.path.isdir(os.path.join(root, folder, crystal_folder)):

            residue_csv = os.path.join(
                root, folder, crystal_folder, "residue_scores.csv"
            )

            if not os.path.isfile(residue_csv):

                os.chdir(os.path.join(root, folder, crystal_folder))
                os.system(
                    "giant.score_model {} {}".format(
                        folders[folder][0], folders[folder][1]
                    )
                )

            if not os.path.isfile(residue_csv):
                continue

            scores_df = pd.read_csv(residue_csv)
            scores_df["crystal"] = crystal_folder
            scores_df["method"] = folders[folder][2]
            scores_df["source"] = folder
            scores_df_list.append(scores_df)

os.chdir("/dls/science/groups/i04-1/elliot-dev/Work/NUDT7A_mass_spec_refinements/")
scores_df = pd.concat(scores_df_list)
scores_df.to_csv("residue_scores.csv")
