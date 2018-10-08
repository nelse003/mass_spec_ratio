import pandas as pd
import os
import csv
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from itertools import chain

#"/hdlocal/home/enelson/mass_spec_ratio/Data/NUDT7-p019-70904_h3d-0_4.CSV"

def get_peak_range(df):

    start= df['Start X'].min()
    end = df['End X'].max()
    range = end - start
    peaks_x_range = (df['End X'] - df['Start X']).sum()

    return range, peaks_x_range

def get_ratios(csv_name, fnames):

    ms_deconv_df = pd.read_csv(csv_name, header=2, index_col=False)
    ms_deconv_df = ms_deconv_df.dropna(axis='columns', how='all')
    ms_deconv_df = ms_deconv_df[ms_deconv_df['Center X'].between(20000,30000, inclusive=True)]

    unlabelled_df = ms_deconv_df[ms_deconv_df['Center X'].between(25100,25350, inclusive=True)]
    labelled_df = ms_deconv_df[ms_deconv_df['Center X'].between(25400,25650, inclusive=True)]


    height_ranked_df = ms_deconv_df.sort_values(by=['Height'],ascending=False).reset_index(drop=True)

    if not unlabelled_df.empty:
        unlabelled_peak = unlabelled_df.loc[unlabelled_df['Height'].idxmax()]['Center X']
        range_unlabelled, peaks_x_range_unlabelled = get_peak_range(unlabelled_df)
        rank_height_unlabelled_df = height_ranked_df[
            height_ranked_df['Center X'].between(25100, 25350, inclusive=True)]
        rank_height_unlabelled = (rank_height_unlabelled_df['Height']).idxmax()


    if not labelled_df.empty:
        labelled_peak = labelled_df.loc[labelled_df['Height'].idxmax()]['Center X']
        range_labelled, peaks_x_range_labelled = get_peak_range(labelled_df)
        rank_height_labelled_df = height_ranked_df[
            height_ranked_df['Center X'].between(25400, 25650, inclusive=True)]

        rank_height_labelled = (rank_height_labelled_df['Height']).idxmax()

    if unlabelled_df.empty and not labelled_df.empty:
        area_ratio = 1
        peak_area_ratio = 1
        height_ratio = 1
        unlabelled_peak = np.nan
        range_unlabelled = np.nan
        peaks_x_range_unlabelled = np.nan
        rank_height_unlabelled = np.nan

    elif not unlabelled_df.empty and labelled_df.empty:
        area_ratio = 0
        peak_area_ratio = 0
        height_ratio = 0
        labelled_peak = np.nan
        range_labelled = np.nan
        peaks_x_range_labelled = np.nan
        rank_height_labelled = np.nan

    elif unlabelled_df.empty and labelled_df.empty:
        return [csv_name]

    elif not unlabelled_df.empty and not labelled_df.empty:
        area_ratio = labelled_df['Area'].sum()/(unlabelled_df['Area'].sum() + labelled_df['Area'].sum())
        peak_area_ratio = labelled_df['Area'].max()/(unlabelled_df['Area'].max() + labelled_df['Area'].max())
        height_ratio = labelled_df['Height'].max()/(unlabelled_df['Height'].max() + labelled_df['Height'].max())

    return [csv_name, area_ratio, height_ratio, peak_area_ratio,
            unlabelled_peak, labelled_peak, range_unlabelled,
            range_labelled, peaks_x_range_unlabelled, peaks_x_range_labelled,
            rank_height_unlabelled, rank_height_labelled ]


def get_intended_ratio(csv_name):

    csv = os.path.basename(csv_name)
    #ratio = float(csv.split('-')[-1][0:3].replace('_','.'))
    ratio = csv.split('-')[-1][0:3].replace('_','.')
    print(ratio)
    return ratio

data_dir = "/hdlocal/home/enelson/mass_spec_ratio/Data/Solution/"
output_file = "/hdlocal/home/enelson/mass_spec_ratio/Output/solution_ratios.csv"

fnames = ['csv_names','area_ratio', 'height_ratio', 'peak_area_ratio',
'unlabelled_peak', 'labelled_peak', 'range_unlabelled',
'range_labelled', 'peaks_x_range_unlabelled', 'peaks_x_range_labelled',
          'rank_height_unlabelled', 'rank_height_labelled','intended_ratio']

with open(output_file,'w') as csvfile:
    writer = csv.writer(csvfile, delimiter=',')
    writer.writerow(fnames)
    file_gen = (f for f in os.listdir(data_dir) if os.path.isfile(os.path.join(data_dir, f)))
    # low_file_gen = (f for f in os.listdir(os.path.join(data_dir,"Low_signal"))
    #             if os.path.isfile(os.path.join(data_dir,"Low_signal", f)))

    for file in file_gen:
        intended_ratio = get_intended_ratio(file)
        ratios = get_ratios(os.path.join(data_dir, file), fnames)
        ratios.append(intended_ratio)
        writer.writerow(ratios)

    # for file in low_file_gen:
    #     intended_ratio = get_intended_ratio(file)
    #     ratios = get_ratios(os.path.join(data_dir,"Low_signal", file), fnames)
    #     ratios.append(intended_ratio)
    #     writer.writerow(ratios)

df = pd.read_csv(output_file, header=0)

plt.scatter(x=df['intended_ratio'], y=df['height_ratio'], label='Largest Peak Height')
# plt.plot(np.unique(df['intended_ratio']),
#          np.poly1d(np.polyfit(df['intended_ratio'],
#                               df['height_ratio'], 1))
#          (np.unique(df['intended_ratio'])))

plt.scatter(x=df['intended_ratio'], y=df['area_ratio'], label='Area of all adduct peaks')
# plt.plot(np.unique(df['intended_ratio']),
#          np.poly1d(np.polyfit(df['intended_ratio'],
#                               df['area_ratio'], 1))
#          (np.unique(df['intended_ratio'])))

plt.scatter(x=df['intended_ratio'], y=df['peak_area_ratio'], label='Area of largest peak')
# plt.plot(np.unique(df['intended_ratio']),
#          np.poly1d(np.polyfit(df['intended_ratio'],
#                               df['peak_area_ratio'], 1))
#          (np.unique(df['intended_ratio'])))
plt.plot([0,1],[0,1], label='x=y')

plt.xlabel('Intended Ratio')
plt.ylabel('Mass Spec Ratio')
plt.title('In Solution')
plt.legend(loc='best')
plt.savefig("/hdlocal/home/enelson/mass_spec_ratio/Output/mass_spec_ratios_solution.png", dpi=300)