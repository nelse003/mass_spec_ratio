import pandas as pd
import os

def filter_csv(dataset_folder, masses):

    df_list = []
    for csv_name in os.listdir(dataset_folder):

        if csv_name.endswith(".csv"):

            ms_df = pd.read_csv(os.path.join(dataset_folder, csv_name), header=1, index_col=False)
            filter_df = ms_df[ms_df["X(Daltons)"].isin(masses)]
            filter_df['file'] = csv_name

            df_list.append(filter_df)


    filtered_df = pd.concat(df_list)

    filtered_df.to_csv(os.path.join(upper_folder, "{}.csv".format(os.path.basename(dataset_folder))))

SMAD1_c001 = [24735, 24815, 24895]
ACVR1A_c076 = [39785, 39865, 39945,40025,40105,40185,40265,40345,40425]
ACVR1Z_r206h = [37334, 37414, 37494, 37574, 37654, 37734, 37814, 37894, 37974]
ACVR1Z_q207e = [38469, 38549, 38629,	38709,38789,38869,38949,39029,39109]
ACVR2A_c046 = [36409, 36489, 36569,	36649, 36729, 36809, 36889,	36969, 37049, 33975,
               34055, 34135, 34215,	34295, 34375, 34455, 34535,	34615]
BMPR2A_c015 = [40114,40194,40274,40354,40434,40514,40594,40674,40754]
BMPR2Z_c002 = [37634,37714,37794,37874,37954,38034,38114,38194,38274]

upper_folder = "/scratch/Ellie/mass_spec/181031_phos_assay_timecourse/"

Q207E_A2 = SMAD1_c001 + ACVR1Z_q207e + ACVR2A_c046
Q207E_B2A = SMAD1_c001 + ACVR1Z_q207e + BMPR2A_c015
Q207E_B2Z = SMAD1_c001 + ACVR1Z_q207e + BMPR2Z_c002

R206H_A2 = SMAD1_c001 + ACVR1Z_r206h + ACVR2A_c046
R206H_B2A = SMAD1_c001 + ACVR1Z_r206h + BMPR2A_c015
R206H_B2Z = SMAD1_c001 + ACVR1Z_r206h + BMPR2Z_c002

WT_A2 = SMAD1_c001 + ACVR1A_c076  + ACVR2A_c046
WT_B2A = SMAD1_c001 + ACVR1A_c076  + BMPR2A_c015
WT_B2Z = SMAD1_c001 + ACVR1A_c076  + BMPR2Z_c002

for folder in os.listdir(upper_folder):
    if "_" in folder:
        dataset_folder = os.path.join(upper_folder,folder)
        filter_csv(dataset_folder,eval(folder))







