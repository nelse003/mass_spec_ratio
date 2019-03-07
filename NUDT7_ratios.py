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

def get_ratios_of_expected_peaks(df,
                                 expected_unlabelled_peaks,
                                 expected_labelled_peaks):

    ratio_dict = {}
    for peak, peak_labelled in zip(expected_unlabelled_peaks,
                                   expected_labelled_peaks):

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

def ratios_from_filenames(df_dict):

    intended_ratio_dict = {}
    before_crystal_ratio_dict = {}
    for key, df in df_dict.items():

        A_D = ['_A','_B','_C','_D']
        E_H = ['_E', '_F', '_G', '_H']

        if 'CI074433' in key:
            if any(s in key for s in A_D):
                intended_ratio_dict[key] = 0
            if any(s in key for s in E_H):
                intended_ratio_dict[key] = 0.10
        elif 'CI074436' in key:
            if any(s in key for s in A_D):
                intended_ratio_dict[key] = 0.20
            if any(s in key for s in E_H):
                intended_ratio_dict[key] = 0.30
        elif 'CI074435' in key:
            if any(s in key for s in A_D):
                intended_ratio_dict[key] = 0.40
            if any(s in key for s in E_H):
                intended_ratio_dict[key] = 0.50
        elif 'CI074434' in key:
            if any(s in key for s in A_D):
                intended_ratio_dict[key] = 0.60
            if any(s in key for s in E_H):
                intended_ratio_dict[key] = 0.70
        elif 'CI074438' in key:
            if any(s in key for s in A_D):
                intended_ratio_dict[key] = 0.80
            if any(s in key for s in E_H):
                intended_ratio_dict[key] = 0.90
        elif 'CI074437' in key:
            if any(s in key for s in A_D):
                intended_ratio_dict[key] = 1.00
            if any(s in key for s in E_H):
                intended_ratio_dict[key] = 0.75
        elif '0L_100U' in key:
            intended_ratio_dict[key] = 0
        elif '10L_90U' in key:
            intended_ratio_dict[key] = 0.10
        elif '20L_80U' in key:
            intended_ratio_dict[key] = 0.20
        elif '30L_70U' in key:
            intended_ratio_dict[key] = 0.30
        elif '40L_60U' in key:
            intended_ratio_dict[key] = 0.40
        elif '50L_50U' in key:
            intended_ratio_dict[key] = 0.50
        elif '60L_40U' in key:
            intended_ratio_dict[key] = 0.60
        elif '70L_30U' in key:
            intended_ratio_dict[key] = 0.70
        elif '80L_20U' in key:
            intended_ratio_dict[key] = 0.80
        elif '90L_10U' in key:
            intended_ratio_dict[key] = 0.90
        elif '100L_0U' in key:
            intended_ratio_dict[key] = 1.00
        else:
            raise ValueError("Key not recognised: {}".format(key))

    return intended_ratio_dict

def remove_blank_datasets(df_dict):

    keys_to_remove = []
    for key, df in df_dict.items():
        if "blank" in key.lower():
            keys_to_remove.append(key)
    for key in set(keys_to_remove):
        del df_dict[key]

    return df_dict

def remove_dataset_by_filename_content(df_dict, key_string):

    keys_to_remove = []
    for key, df in df_dict.items():
         if key_string.lower() in key.lower():
             keys_to_remove.append(key)
    for key in set(keys_to_remove):
        del df_dict[key]

    return df_dict


def check_signal_to_noise(df_dict):

    keys_to_remove = []
    peak_ratio = []
    for key, df in df_dict.items():
        try:
            peak_ratio.append(peak_height_outside_of_interest(df))
        except ValueError:
            keys_to_remove.append(key)
    for key in set(keys_to_remove):
        del df_dict[key]

    return df_dict, peak_ratio

def string_contains(check_str, match):
    if match in check_str:
        return True
    else:
        return False

def pre_crystal_plot(df):

    x = df['intended_ratio']

    y_h = df['weighted_height_ratio']
    y_a = df['weighted_area_ratio']

    fig, ax = plt.subplots()

    ax.scatter(x, y_h, color='red', label='Height Ratio')
    ax.scatter(x, y_a, color='blue', label='Area Ratio')
    ax.legend()

    plt.ylabel('Measured Ratio of Labelled Species')
    plt.xlabel('Intended Ratio of Labelled Species')
    plt.title('Pre Crystallisation Calibration Curve')
    plt.savefig('Output/pre_crystal_ratio.png')

def post_crystal_plot(df):

    df = df.convert_objects(convert_numeric=True)

    x = np.array(df['intended_ratio'].unique())
    x.sort()

    y_h_mean = df.groupby('intended_ratio')['weighted_height_ratio'].mean()
    y_a_mean = df.groupby('intended_ratio')['weighted_area_ratio'].mean()
    y_a_std = df.groupby('intended_ratio')['weighted_area_ratio'].std()
    y_h_std = df.groupby('intended_ratio')['weighted_height_ratio'].std()

    print(x)
    print(y_h_mean)

    fig, ax = plt.subplots()

    ax.errorbar(x=x, y=y_h_mean, yerr=y_h_std, color='red',fmt='o', label='Height Ratio (std)')
    ax.errorbar(x=x, y=y_a_mean, yerr=y_a_std, color='blue',fmt='o', label='Area Ratio (std)')
    ax.legend()

    plt.ylabel('Measured Ratio of Labelled Species')
    plt.xlabel('Intended Ratio of Labelled Species')
    plt.title('Post Crystallisation Calibration Curve')
    plt.savefig('Output/post_crystal_ratio.png')

def post_crystal_calibrated_plot(df):

    df = df.convert_objects(convert_numeric=True)

    x = np.array(df['calibrated_ratio'].unique())
    x.sort()

    y_h_mean = df.groupby('calibrated_ratio')['weighted_height_ratio'].mean()
    y_a_mean = df.groupby('calibrated_ratio')['weighted_area_ratio'].mean()
    y_a_std = df.groupby('calibrated_ratio')['weighted_area_ratio'].std()
    y_h_std = df.groupby('calibrated_ratio')['weighted_height_ratio'].std()

    fig, ax = plt.subplots()

    ax.errorbar(x=x, y=y_h_mean, yerr=y_h_std, color='red',fmt='o', label='Height Ratio (std)')
    ax.errorbar(x=x, y=y_a_mean, yerr=y_a_std, color='blue',fmt='o', label='Area Ratio (std)')
    ax.legend()

    plt.ylabel('Measured Ratio of Labelled Species')
    plt.xlabel('Calibrated Ratio of Labelled Species')
    plt.title('Post Crystallisation with Calibrated Ratio')
    plt.savefig('Output/post_crystal_ratio.png')

if __name__ == "__main__":

    """
    Notes
    ------
    conda create --name mass_spec_ratio numpy scipy matplotlib pandas notebook python=3.6
    """

    #data_dir = "/hdlocal/enelson/mass_spec_ratio/NUDT7_Data"
    #data_dir = "/home/nelse003/PycharmProjects/mass_spec_ratio/NUDT7_Data"
    data_dir = "/dls/science/groups/i04-1/elliot-dev/mass_spec_ratio/NUDT7_Data"

    # Parse all csv in folder to separate into individual deconvolutions by headers
    df_dict = {}
    for csv in os.listdir(data_dir):
        csv = os.path.join(data_dir, csv)
        df_dict.update(read_grouped_csv(csv))

    # Remove blank datasets
    df_dict = remove_dataset_by_filename_content(df_dict, key_string="blank")
    # Remove un-needed dataset
    df_dict = remove_dataset_by_filename_content(df_dict,
                                                 key_string="NUDT7A_p026_NU0000308a"
                                                            "_post_gel_filtration")

    # remove datasets that do not pass signal threshold
    # currently 1000 counts in interest area 22000-26000
    # TODO Move parameters up to non hard coded region
    df_dict, peak_ratios = check_signal_to_noise(df_dict)

    # From Exploratory analysis these pairs of peaks are present
    # showing that degradation into three sets of paired peaks
    expected_unlabelled_peaks = np.array([25125, 24218, 23525, 22786])
    expected_labelled_peaks = expected_unlabelled_peaks + 354

    # Translate keys (filenames) into ratios
    # Split across two dicitionaries to parse
    # intended ratio and expected ratio
    intended_ratio_dict = ratios_from_filenames(df_dict)

    # Process all deconvolutions to ratios of peaks
    ratio_df_list = []
    for key, df in df_dict.items():
        ratio_df = get_ratios_of_expected_peaks(df,
                                                expected_unlabelled_peaks,
                                                expected_labelled_peaks)

        # Get weights of contributing peaks based on peak height
        height_weight = (ratio_df['peak_height'] + ratio_df['labelled_peak_height']) \
                        / (ratio_df['labelled_peak_height'].sum() + \
                           ratio_df['peak_height'].sum())

        # Get weights of contributing peaks based on peak area
        area_weight = (ratio_df['peak_area'] + ratio_df['labelled_peak_area']) \
                      / (ratio_df['labelled_peak_area'].sum() + \
                         ratio_df['peak_area'].sum())

        # Get a single weighted ratio using peak heights
        weighted_height_ratio = height_weight * ratio_df['height_ratio']
        weighted_height_ratio = weighted_height_ratio.sum()

        # Get a single weighted ratio using peak areas
        weighted_area_ratio = area_weight * ratio_df['area_ratio']
        weighted_area_ratio = weighted_area_ratio.sum()

        # Store weights and ratios, and key in ratio dataframe
        ratio_df['height_weights'] = height_weight
        ratio_df['area_weights'] = area_weight
        ratio_df['weighted_height_ratio'] = weighted_height_ratio
        ratio_df['weighted_area_ratio'] = weighted_area_ratio
        ratio_df['intended_ratio'] = intended_ratio_dict[key]
        ratio_df['key'] = key

        # For concatenating results
        ratio_df_list.append(ratio_df)

    ratio_df = pd.concat(ratio_df_list)

    # Split into pre crystallisation and post crystallisation
    ratio_df['pre_crystal'] = ratio_df['key'].apply(string_contains, match='L_')
    pre_crystal_df = ratio_df[ratio_df['pre_crystal'] == True]
    post_crystal_df = ratio_df[ratio_df['pre_crystal'] == False]

    intended_to_weighted = dict(zip(pre_crystal_df['intended_ratio'],
                                    pre_crystal_df['weighted_height_ratio']))

    pre_crystal_plot(pre_crystal_df)
    post_crystal_plot(post_crystal_df)

    # Add calibration from pre crystal to post crystal
    df = pre_crystal_df[['intended_ratio', 'weighted_height_ratio']]
    df = df.drop_duplicates()
    df = df.sort_values(by=['intended_ratio'])
    intended_ratio = df['intended_ratio']
    expected_ratio = df['weighted_height_ratio']
    interp_val = np.interp(0.75, intended_ratio, expected_ratio)

    df1 = pd.DataFrame({'intended_ratio': [0.75, 0],
                        'weighted_height_ratio': [interp_val, 0]})
    df = df.append(df1)
    df = df.sort_values(by=['intended_ratio'])
    df = df.rename(columns={'weighted_height_ratio': 'calibrated_ratio'})
    post_crystal_df = pd.merge(post_crystal_df, df, on='intended_ratio')

    post_crystal_calibrated_plot(post_crystal_df)