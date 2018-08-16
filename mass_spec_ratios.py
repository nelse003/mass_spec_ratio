import pandas as pd
import os
import csv
import numpy as np

#"/hdlocal/home/enelson/mass_spec_ratio/Data/NUDT7-p019-70904_h3d-0_4.CSV"

def get_ratios(csv_name, fnames):

    print(csv_name)

    ms_deconv_df = pd.read_csv(csv_name,header=2, index_col=False)
    ms_deconv_df = ms_deconv_df.dropna(axis='columns', how='all')
    unlabelled_df = ms_deconv_df[ms_deconv_df['Center X'].between(25100,25350, inclusive=True)]
    labelled_df = ms_deconv_df[ms_deconv_df['Center X'].between(25400,25650, inclusive=True)]

    if not unlabelled_df.empty:
        unlabelled_peak = unlabelled_df.loc[unlabelled_df['Height'].idxmax()]['Center X']

        # Find range of peaks
        start_range_unlabelled = unlabelled_df['Start X'].min()
        end_range_unlabelled = unlabelled_df['End X'].max()
        range_unlabelled = end_range_unlabelled - start_range_unlabelled
        peaks_x_range_unlabelled = (unlabelled_df['End X'] - unlabelled_df['Start X']).sum()

    if not labelled_df.empty:
        labelled_peak = labelled_df.loc[labelled_df['Height'].idxmax()]['Center X']

        start_range_labelled = labelled_df['Start X'].min()
        end_range_labelled = labelled_df['End X'].max()
        range_labelled = end_range_labelled - start_range_labelled
        peaks_x_range_labelled = (labelled_df['End X'] - labelled_df['Start X']).sum()

    if unlabelled_df.empty and not labelled_df.empty:
        area_ratio = 1
        peak_area_ratio = 1
        height_ratio = 1
        unlabelled_peak = np.nan
        range_unlabelled = np.nan
        peaks_x_range_unlabelled = np.nan

    elif not unlabelled_df.empty and labelled_df.empty:
        area_ratio = 0
        peak_area_ratio = 0
        height_ratio = 0
        labelled_peak = np.nan
        range_labelled = np.nan
        peaks_x_range_labelled = np.nan

    elif unlabelled_df.empty and labelled_df.empty:
        return [np.nan * len(fnames)]

    elif not unlabelled_df.empty and not labelled_df.empty:
        area_ratio = labelled_df['Area'].sum()/unlabelled_df['Area'].sum()
        peak_area_ratio = labelled_df['Area'].max()/unlabelled_df['Area'].max()
        height_ratio = labelled_df['Height'].max()/unlabelled_df['Height'].max()

    return [area_ratio, height_ratio, peak_area_ratio,
            unlabelled_peak, labelled_peak, range_unlabelled,
            range_labelled, peaks_x_range_unlabelled, peaks_x_range_labelled]

fnames = ['area_ratio', 'height_ratio', 'peak_area_ratio',
'unlabelled_peak', 'labelled_peak', 'range_unlabelled',
'range_labelled', 'peaks_x_range_unlabelled', 'peaks_x_range_labelled']

data_dir = "/hdlocal/home/enelson/mass_spec_ratio/Data/"

with open("/hdlocal/home/enelson/mass_spec_ratio/Output/ratios.csv",'w') as csvfile:
    writer = csv.writer(csvfile, delimiter=',')
    writer.writerow(fnames)
    for file in os.listdir(data_dir):
        ratios = get_ratios(os.path.join(data_dir,file), fnames)
        writer.writerow(ratios)