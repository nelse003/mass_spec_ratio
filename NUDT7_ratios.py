import os
import pandas as pd
import numpy as np
from scipy.signal import find_peaks
import re
from io import StringIO

import matplotlib.pyplot as plt
plt.rcParams['figure.figsize'] = [10, 10]

pd.set_option("display.max_columns", 10)

def grouped(iterable, n):
    """
    Chunking iterable into sizes n

    Notes
    ----------
    Python3 uses zip, python 2 would require izip from itertools

    Examples
    ----------
    s -> (s0,s1,s2,...sn-1), (sn,sn+1,sn+2,...s2n-1), (s2n,s2n+1,s2n+2,...s3n-1), ...
    """
    return zip(*[iter(iterable)]*n)


def read_grouped_csv(csv):

    """
    Open a csv split by headers preceeed by #

    Ouputs a dictionary of dataframes, prefixed by the header.

    :param csv:
    :return:
    """
    date = os.path.basename(csv).split('_')[0]

    with open(csv) as f:
        data = f.read()

    header_re = re.compile("#*")
    split_data = header_re.split(data)
    #if not 'Point' in split_data[0]:
    split_data.pop(0)

    df_dict = {}
    for header, data in grouped(split_data, 2):
        df_dict[date + " : " + header] = pd.read_csv(StringIO(data))

    return df_dict


def peak_height_outside_of_interest(df, min_interest=22000, max_interest=26000):
    """ Find ratio between largest peak inside area of interest, and outside """

    peaks, _ = find_peaks(df['Y(Counts)'], height=1000, distance=100)

    peaks_df = df.loc[peaks]

    interest_peaks_df = peaks_df[(peaks_df['X(Daltons)'] >= min_interest)
                                 & (peaks_df['X(Daltons)'] <= max_interest)]
    uninterest_peaks_df = peaks_df[(peaks_df['X(Daltons)'] < min_interest)
                                   | (peaks_df['X(Daltons)'] > max_interest)]

    if interest_peaks_df.empty:
       raise ValueError("Signal is below cutoff threshold."
                        " Data is too weak to be analysed")
    elif uninterest_peaks_df.empty:
       raise ValueError("Signal is below cutoff threshold."
                        " Data is too weak to be analysed")


    max_uninterest_peak = uninterest_peaks_df.nlargest(1, 'Y(Counts)')
    max_interest_peak = interest_peaks_df.nlargest(1, 'Y(Counts)')

    ratio_of_peak_outside_area = max_uninterest_peak['Y(Counts)'].values / (
                max_interest_peak['Y(Counts)'].values +
                max_uninterest_peak['Y(Counts)'].values)

    return float(ratio_of_peak_outside_area)

def peak_nearby(df, peak, delta=10, left_range=50, right_range=100):

    peak_df = df[(df['X(Daltons)'] >= peak - delta)
                 & (df['X(Daltons)'] <= peak + delta)]
    peak_loc = peak_df.loc[peak_df['Y(Counts)'].idxmax()]['X(Daltons)']
    peak_height = peak_df['Y(Counts)'].max()
    peak_area_df = df[(df['X(Daltons)'] >= peak - left_range)
                    & (df['X(Daltons)'] <= peak + right_range)]
    peak_area = peak_area_df['Y(Counts)'].sum()

    return peak_loc, peak_height, peak_area

def get_ratios_of_expected_peaks(df,expected_unlabelled_peaks,expected_labelled_peaks):

    ratio_dict = {}
    for peak, peak_labelled in zip(expected_unlabelled_peaks, expected_labelled_peaks):

        peak_loc, peak_height, peak_area = peak_nearby(df, peak,
                                                       delta=10,
                                                       left_range=50,
                                                       right_range=100)
        labelled_peak_loc, \
        labelled_peak_height, \
        labelled_peak_area = peak_nearby(df,
                                         peak_labelled,
                                         delta=10,
                                         left_range=50, right_range=100)

        height_ratio = labelled_peak_height/(labelled_peak_height + peak_height)
        area_ratio = labelled_peak_area/(labelled_peak_area + peak_area)
        ratio_dict[peak] = (height_ratio,
                            area_ratio,
                            labelled_peak_loc,
                            labelled_peak_height,
                            labelled_peak_area,
                            peak_loc,
                            peak_height,
                            peak_area)

    ratio_df = pd.DataFrame.from_dict(ratio_dict,
                                      orient='index',
                                      columns=['height_ratio',
                                               'area_ratio',
                                               'labelled_peak_loc',
                                               'labelled_peak_height',
                                               'labelled_peak_area',
                                               'peak_loc',
                                               'peak_height',
                                               'peak_area'
                                               ])

    return ratio_df

if __name__ == "__main__":

    # data_dir = "/hdlocal/enelson/mass_spec_ratio/NUDT7_Data"
    data_dir = "/home/nelse003/PycharmProjects/mass_spec_ratio/NUDT7_Data"

    df_dict = {}
    for csv in os.listdir(data_dir):
        csv = os.path.join(data_dir, csv)
        df_dict.update(read_grouped_csv(csv))

    keys_to_remove = []
    peak_ratio  = []
    for key, df in df_dict.items():
        if "blank" in key:
            keys_to_remove.append(key)

    for key in set(keys_to_remove):
        del df_dict[key]

    keys_to_remove = []
    for key, df in df_dict.items():
        try:
            peak_ratio.append(peak_height_outside_of_interest(df))
        except ValueError:
            keys_to_remove.append(key)

    expected_unlabelled_peaks = np.array([25125, 24218, 23525, 22786])
    expected_labelled_peaks = expected_unlabelled_peaks + 354

    key_ratio_dict = {}
    for key, df in df_dict.items():

        A_D = ['_A','_B','_C','_D']
        E_H = ['_E', '_F', '_G', '_H']

        if 'CI074433' in key:
            if any(s in key for s in A_D):
                key_ratio_dict[key] = 0
            if any(s in key for s in E_H):
                key_ratio_dict[key] = 0.10
        elif 'CI074436' in key:
            if any(s in key for s in A_D):
                key_ratio_dict[key] = 0.20
            if any(s in key for s in E_H):
                key_ratio_dict[key] = 0.30
        elif 'CI074435' in key:
            if any(s in key for s in A_D):
                key_ratio_dict[key] = 0.40
            if any(s in key for s in E_H):
                key_ratio_dict[key] = 0.50
        elif 'CI074434' in key:
            if any(s in key for s in A_D):
                key_ratio_dict[key] = 0.60
            if any(s in key for s in E_H):
                key_ratio_dict[key] = 0.70
        elif 'CI074438' in key:
            if any(s in key for s in A_D):
                key_ratio_dict[key] = 0.80
            if any(s in key for s in E_H):
                key_ratio_dict[key] = 0.90
        elif 'CI074437' in key:
            if any(s in key for s in A_D):
                key_ratio_dict[key] = 1.00
            if any(s in key for s in E_H):
                key_ratio_dict[key] = 0.75
        else:
            print("Fuck: {}".format(key))

    for key, df in df_dict.items():
        ratio_df = get_ratios_of_expected_peaks(df,
                                                expected_unlabelled_peaks,
                                                expected_labelled_peaks)

        height_weights = (ratio_df['peak_height'] + ratio_df['labelled_peak_height']) \
              /(ratio_df['labelled_peak_height'].sum() + \
                ratio_df['peak_height'].sum())

        area_weight = (ratio_df['peak_area'] + ratio_df['labelled_peak_area']) \
              /(ratio_df['labelled_peak_area'].sum() + \
                ratio_df['peak_area'].sum())

        weighted_height_ratio = height_weights * ratio_df['height_ratio']
        weighted_height_ratio = weighted_height_ratio.sum()
        weighted_area_ratio = area_weight * ratio_df['area_ratio']
        weighted_area_ratio = weighted_area_ratio.sum()

        print(key_ratio_dict[key], weighted_height_ratio, weighted_area_ratio)



