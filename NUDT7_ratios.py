import os
import pandas as pd
import numpy as np
from scipy.signal import find_peaks
import re
from io import StringIO

import matplotlib.pyplot as plt
import statsmodels.api as sm
import argparse

plt.rcParams["figure.figsize"] = [10, 10]

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
    return zip(*[iter(iterable)] * n)


def read_grouped_csv(csv):

    """
    Open a csv split by headers preceeed by #

    Ouputs a dictionary of dataframes, prefixed by the header.

    Parameters
    ----------
    csv: str
        path to csv file with multiple mass spectroscopy data

    Returns
    -------
    df_dict: dict
        dictionary of pandas.DataFrames
        each containing mass spectroscopy data:
        Y(Counts) vs X(Daltons) data
    """
    date = os.path.basename(csv).split("_")[0]

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

    """
    Find ratio between largest peak inside area of interest, and outside

     Parameters
     -----------
     df: pandas.DataFrame
        DataFrame containing mass spectroscopy data:
         Y(Counts) vs X(Daltons) data

     min_interest: int or float
        Minima of range to be explored (Daltons)

     max_interest:
        Maxima of range to be explored (Daltons)

     Returns
     --------
     ratio_of_peak_outside_area :float
        ratio between largest peak inside area of interest,
        and outside area of interest.
     """

    # Get peaks
    peaks, _ = find_peaks(df["Y(Counts)"], height=1000, distance=100)

    # Turn peak information into dataframe
    peaks_df = df.loc[peaks]

    # Select peaks in the interesting region
    interest_peaks_df = peaks_df[
        (peaks_df["X(Daltons)"] >= min_interest)
        & (peaks_df["X(Daltons)"] <= max_interest)
    ]

    # Select peaks not in the interesting region
    uninterest_peaks_df = peaks_df[
        (peaks_df["X(Daltons)"] < min_interest)
        | (peaks_df["X(Daltons)"] > max_interest)
    ]

    # If no peaks are found then we assume the signal to noise ratio is too low
    if interest_peaks_df.empty:
        raise ValueError(
            "Signal is below cutoff threshold. Data is too weak to be analysed"
        )

    # Get the index of largest peak values
    max_uninterest_peak = uninterest_peaks_df.nlargest(1, "Y(Counts)")
    max_interest_peak = interest_peaks_df.nlargest(1, "Y(Counts)")

    # Get ratio between the peak in the intersting region
    # the un-interesting region
    ratio_of_peak_outside_area = max_uninterest_peak["Y(Counts)"].values / (
        max_interest_peak["Y(Counts)"].values + max_uninterest_peak["Y(Counts)"].values
    )

    return float(ratio_of_peak_outside_area)


def get_peak_area(df, peak, delta=10, left_range=50, right_range=100):
    """
    Find peak area

    Check for peak within delta of supplied peak location.
    Get area of peak within left_range and right_range

    Parameters
    ----------
    df: pandas.DataFrame
        DataFrame containing mass spectroscopy data:
         Y(Counts) vs X(Daltons) data

    peak: float
        specified peak location

    delta: float
        distance from supplied peak to check
        for a peak within the peak DataFrame

    left_range: float
        distance from found peak to get area leftwards

    right_range: float
        distance from found peak to get area rightwards

    Returns
    -------
    peak_loc:

    peak_height:

    peak_area:

    """
    # Trim the dataframe to data near peak location (within delta)
    peak_df = df[
        (df["X(Daltons)"] >= peak - delta) & (df["X(Daltons)"] <= peak + delta)
    ]

    # Get the index of the peak in the peak dataframe
    peak_loc = peak_df.loc[peak_df["Y(Counts)"].idxmax()]["X(Daltons)"]

    # Get the peak height
    peak_height = peak_df["Y(Counts)"].max()

    # Get the
    peak_area_df = df[
        (df["X(Daltons)"] >= peak - left_range)
        & (df["X(Daltons)"] <= peak + right_range)
    ]

    peak_area = peak_area_df["Y(Counts)"].sum()

    return peak_loc, peak_height, peak_area


def get_ratios_of_expected_peaks(
    df, expected_unlabelled_peaks, expected_labelled_peaks
):

    ratio_dict = {}
    for peak, peak_labelled in zip(expected_unlabelled_peaks, expected_labelled_peaks):

        peak_loc, peak_height, peak_area = get_peak_area(
            df, peak, delta=10, left_range=50, right_range=100
        )
        labelled_peak_loc, labelled_peak_height, labelled_peak_area = get_peak_area(
            df, peak_labelled, delta=10, left_range=50, right_range=100
        )

        height_ratio = labelled_peak_height / (labelled_peak_height + peak_height)
        area_ratio = labelled_peak_area / (labelled_peak_area + peak_area)
        ratio_dict[peak] = (
            height_ratio,
            area_ratio,
            labelled_peak_loc,
            labelled_peak_height,
            labelled_peak_area,
            peak_loc,
            peak_height,
            peak_area,
        )

    ratio_df = pd.DataFrame.from_dict(
        ratio_dict,
        orient="index",
        columns=[
            "height_ratio",
            "area_ratio",
            "labelled_peak_loc",
            "labelled_peak_height",
            "labelled_peak_area",
            "peak_loc",
            "peak_height",
            "peak_area",
        ],
    )

    return ratio_df


def ratios_from_csv(csv, df_dict):

    """
    Get ratios and experiment type from

    Parameters
    ----------
    csv: str
        path to csv file defining the ratio and experiment type
    df_dict: dict
        dictionary of pandas.DataFrames
        each containing mass spectroscopy data:
        Y(Counts) vs X(Daltons) data

    Returns
    -------

    """

    if os.path.exists(csv):
        ratio_df = pd.read_csv(csv)
    else:
        ratio_df = pd.DataFrame(columns=["f_name", "intended_ratio", "experiment"])

    for f_name, df in df_dict.items():
        print(f_name)

        if f_name in ratio_df.f_name.values:
            continue
        else:
            ratio_df = ratio_df.append(
                {"f_name": f_name, "intended_ratio": np.nan, "experiment": "unknown"},
                ignore_index=True,
            )

    ratio_df.to_csv(csv, index=None, header=True)
    return ratio_df


def ratios_from_filenames(df_dict):
    """
    Match the filename to the expected ratio

    Parameters
    ----------
    df_dict: dict
        dictionary of pandas.DataFrames
        each containing mass spectroscopy data:
        Y(Counts) vs X(Daltons) data

    Returns
    -------
    intended_ratio: dict
        dictionary of intended ratios

    df_dict: dict
        dictionary of pandas.DataFrames
        each containing mass spectroscopy data:
        Y(Counts) vs X(Daltons) data
    """

    intended_ratio = {}

    key_to_remove = []
    for key, df in df_dict.items():

        A_D = ["_A", "_B", "_C", "_D"]
        E_H = ["_E", "_F", "_G", "_H"]

        post_diffract_df = pd.read_csv("post_diffraction_ratios.csv")

        if "NUDT7A_Post_diffraction_crystal_mass_spectra_apr_14_16_2019" in key:
            well = key.split("_")[-1]

            intended_ratio[key] = float(
                post_diffract_df.loc[post_diffract_df["Solubilisation Well"] == well][
                    "Ratio"
                ].values[0]
            )

            if np.isnan(intended_ratio[key]):
                key_to_remove.append(key)

        elif "180425" in key:
            ratio_str = key.split("1day")[1]
            if ratio_str == "":
                if "30min" in key:
                    intended_ratio[key] = 1.0
                else:
                    intended_ratio[key] = 0.0
            elif "No_cov" in key:
                intended_ratio[key] = 0.0
            else:
                print(ratio_str)
                print(
                    ratio_str.split("no_cov")[0]
                    .lstrip("_")
                    .rstrip("_")
                    .replace("_", ".")
                )
                intended_ratio[key] = (
                    float(
                        ratio_str.split("no_cov")[0]
                        .lstrip("_")
                        .rstrip("_")
                        .replace("_", ".")
                    )
                    / 10
                )

        elif "CI074433" in key:
            if any(s in key for s in A_D):
                intended_ratio[key] = 0
            if any(s in key for s in E_H):
                intended_ratio[key] = 0.10

        elif "CI074436" in key:
            if any(s in key for s in A_D):
                intended_ratio[key] = 0.20
            if any(s in key for s in E_H):
                intended_ratio[key] = 0.30

        elif "CI074435" in key:0.8
            if any(s in key for s in A_D):
                intended_ratio[key] = 0.40
            if any(s in key for s in E_H):
                intended_ratio[key] = 0.50

        elif "CI074434" in key:
            if any(s in key for s in A_D):
                intended_ratio[key] = 0.60
            if any(s in key for s in E_H):
                intended_ratio[key] = 0.70

        elif "CI074438" in key:
            if any(s in key for s in A_D):
                intended_ratio[key] = 0.80
            if any(s in key for s in E_H):
                intended_ratio[key] = 0.90

        elif "CI074437" in key:
            if any(s in key for s in A_D):
                intended_ratio[key] = 1.00
            if any(s in key for s in E_H):
                intended_ratio[key] = 0.75

        elif "0L_100U" in key:
            intended_ratio[key] = 0

        elif "10L_90U" in key:
            intended_ratio[key] = 0.10

        elif "20L_80U" in key:
            intended_ratio[key] = 0.20

        elif "30L_70U" in key:
            intended_ratio[key] = 0.30

        elif "40L_60U" in key:
            intended_ratio[key] = 0.40

        elif "50L_50U" in key:
            intended_ratio[key] = 0.50

        elif "60L_40U" in key:
            intended_ratio[key] = 0.60

        elif "70L_30U" in key:
            intended_ratio[key] = 0.70

        elif "80L_20U" in key:
            intended_ratio[key] = 0.80

        elif "90L_10U" in key:
            intended_ratio[key] = 0.90

        elif "100L_0U" in key:
            intended_ratio[key] = 1.00

        elif "0_0" in key:
            intended_ratio[key] = 0.00

        elif "0_1" in key:
            intended_ratio[key] = 0.076

        elif "0_2" in key:
            intended_ratio[key] = 0.145

        elif "0_3" in key:
            intended_ratio[key] = 0.226

        elif "0_4" in key:
            intended_ratio[key] = 0.313

        elif "0_5" in key:
            intended_ratio[key] = 0.406

        elif "0_6" in key:
            intended_ratio[key] = 0.506

        elif "0_7" in key:
            intended_ratio[key] = 0.614

        elif "0_8" in key:
            intended_ratio[key] = 0.732

        elif "0_9" in key:
            intended_ratio[key] = 0.86

        elif "0_95" in key:
            intended_ratio[key] = 0.92

        elif "1_0" in key:
            intended_ratio[key] = 1.00

        else:
            raise ValueError("Key not recognised: {}".format(key))

    for key in set(key_to_remove):
        del df_dict[key]
        del intended_ratio[key]

    return intended_ratio, df_dict


def remove_dataset_by_filename_content(df_dict, key_string):
    """
    Remove datasets that have matching text in their key_string

    Parameters
    ----------
    df_dict: dict
        dictionary of pandas.DataFrames
        each containing mass spectroscopy data:
        Y(Counts) vs X(Daltons) data

    key_string: str
        string describing a dataset to be removed

    Returns
    -------
    df_dict: dict
        dictionary of pandas.DataFrames
        each containing mass spectroscopy data:
        Y(Counts) vs X(Daltons) data
    """

    # Iterate over all datasets
    # to check for datasets that match supplied string
    keys_to_remove = []
    for key, df in df_dict.items():
        if key_string.lower() in key.lower():
            keys_to_remove.append(key)

    # Remove datasets which match supplied key_string
    for key in set(keys_to_remove):
        del df_dict[key]

    return df_dict


def remove_datasets_with_low_signal_to_noise(df_dict):
    """
    Remove datasets with low signal to noise

    Parameters
    ----------
    df_dict: dict
        dictionary of pandas.DataFrames
        each containing mass spectroscopy data:
        Y(Counts) vs X(Daltons) data

    Returns
    -------
    df_dict: dict
        dictionary of pandas.DataFrames
        each containing mass spectroscopy data:
        Y(Counts) vs X(Daltons) data
    """

    # Iterate over all datasets
    # to check for datasets with low signal to noise
    keys_to_remove = []
    for key, df in df_dict.items():
        try:
            peak_height_outside_of_interest(df)
        except ValueError:
            keys_to_remove.append(key)

    # Remove any low signal to noise datasets from df_dict
    for key in set(keys_to_remove):
        del df_dict[key]

    return df_dict


def string_contains(check_str, match):
    """
    Check that a string contains another string

    Parameters
    ----------
    check_str: str
        string to checked
    match: str or list
        str to check for

    Returns
    -------
    bool
        True if match is in check_str
        False if match is not in check_str
    """
    if any(m in check_str for m in match):
        return True
    else:
        return False


def marker_match(row, match, marker):
    """
    Check that a string contains another string

    Parameters
    ----------
    check_str: str
        string to checked
    match: str or list
        str to check for

    Returns
    -------
    marker: str
    """
    print(row["key"])
    if any(m in row["key"] for m in match):
        return marker
    else:
        return row["marker"]


def pre_crystal_plot(df):

    df["marker"] = "o"

    date_m = {
        "190402": "x",
        "190220": "*",
        "180425": "+",
        "190225": "v",
        "190506": "s",
        "imp": "1",
    }

    for date, marker in date_m.items():
        df["marker"] = df.apply(marker_match, match=[date], marker=marker, axis=1)

    fig, ax = plt.subplots()

    summary_df = pd.DataFrame()
    for date, marker in date_m.items():

        # To plot just recent stuff
        if date == "190506":
            pass
        elif date == "imp":
            pass
        else:
            continue

        df_plot = df[df["marker"] == marker]
        x = df_plot["intended_ratio"]
        y_h = df_plot["weighted_height_ratio"]
        y_a = df_plot["weighted_area_ratio"]
        ax.scatter(
            x, y_h, color="red", label="Height Ratio: {}".format(date), marker=marker
        )
        ax.scatter(
            x, y_a, color="blue", label="Area Ratio: {}".format(date), marker=marker
        )
        if summary_df.empty:
            summary_df = df_plot
        else:
            summary_df = pd.concat([summary_df, df_plot])

    fit_area = sm.OLS(
        summary_df["weighted_area_ratio"], sm.add_constant(summary_df["intended_ratio"])
    ).fit()

    print(fit_area.summary())

    print(fit_area.params[0])

    plt.plot(
        np.linspace(0, 1, 100),
        np.linspace(0, 1, 100) * fit_area.params[1] + fit_area.params[0],
        label="y={}x+{}".format(fit_area.params[1], fit_area.params[0]),
    )

    plt.legend()
    plt.ylabel("Measured Ratio of Labelled Species")
    plt.xlabel("Intended Ratio of Labelled Species")
    plt.title("Pre Crystallisation Calibration Curve")
    plt.savefig("Output/pre_crystal_ratio.png")


def post_crystal_plot(df):

    df = df.convert_objects(convert_numeric=True)

    x = np.array(df["intended_ratio"].unique())
    x.sort()

    y_h_mean = df.groupby("intended_ratio")["weighted_height_ratio"].mean()
    y_a_mean = df.groupby("intended_ratio")["weighted_area_ratio"].mean()
    y_a_std = df.groupby("intended_ratio")["weighted_area_ratio"].std()
    y_h_std = df.groupby("intended_ratio")["weighted_height_ratio"].std()

    fig, ax = plt.subplots()

    ax.errorbar(
        x=x, y=y_h_mean, yerr=y_h_std, color="red", fmt="o", label="Height Ratio (std)"
    )
    ax.errorbar(
        x=x, y=y_a_mean, yerr=y_a_std, color="blue", fmt="o", label="Area Ratio (std)"
    )

    ax.legend()

    plt.ylabel("Measured Ratio of Labelled Species")
    plt.xlabel("Intended Ratio of Labelled Species")
    plt.title("Post Crystallisation Calibration Curve")
    plt.savefig("Output/post_crystal_ratio.png")


def calibrated_plot(df, ax, marker, label):

    x = np.array(df["calibrated_ratio"].unique())
    x.sort()

    y_h_mean = df.groupby("calibrated_ratio")["weighted_height_ratio"].mean()
    y_a_mean = df.groupby("calibrated_ratio")["weighted_area_ratio"].mean()
    y_a_std = df.groupby("calibrated_ratio")["weighted_area_ratio"].std()
    y_h_std = df.groupby("calibrated_ratio")["weighted_height_ratio"].std()

    ax.errorbar(
        x=x,
        y=y_h_mean,
        yerr=y_h_std,
        color="red",
        fmt=marker,
        label="Height Ratio (std): {}".format(label),
    )
    ax.errorbar(
        x=x,
        y=y_a_mean,
        yerr=y_a_std,
        color="blue",
        fmt=marker,
        label="Area Ratio (std): {}".format(label),
    )

    return ax


def post_crystal_calibrated_plot(df, diffract_df):

    df = df.convert_objects(convert_numeric=True)

    fig, ax = plt.subplots()
    ax = calibrated_plot(df, ax, marker="o", label=" ")
    ax = calibrated_plot(diffract_df, ax, marker="*", label="post diffraction")

    ax.legend()

    plt.ylabel("Measured Ratio of Labelled Species")
    plt.xlabel("Calibrated Ratio of Labelled Species")
    plt.title("Post Crystallisation with Calibrated Ratio")
    plt.savefig("Output/post_crystal_ratio.png")


if __name__ == "__main__":

    """
    Parse csv containign the deconvolutions for mass spectroscopy 
    data into plots to look at pre-crystallisation calibration
    
    Notes
    ------
    conda create --name mass_spec_ratio numpy scipy\
     matplotlib pandas statsmodels notebook python=3.6
    """

    # Add ability to parse command line arguments
    parser = argparse.ArgumentParser(description="Plot NUDT7 mass spectroscopy ratios")
    # Add path argument
    parser.add_argument("--path", action="store", dest="path")
    # resolve the parser
    args = parser.parse_args()

    # Get file paths from input path
    data_dir = os.path.join(args.path, "NUDT7_Data")
    output_dir = os.path.join(args.path, "Output")

    # Create output directory if missing
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # Parse csv and folders of csv.
    df_dict = {}

    # If folder, expect these csv to represent one dataset
    for f in os.listdir(data_dir):
        print(f)

        if os.path.isdir(os.path.join(data_dir, f)):
            for csv in os.listdir(os.path.join(data_dir, f)):
                print(csv)
                df_dict.update(
                    {
                        "{}_{}".format(f, csv.rstrip(".CSV")): pd.read_csv(
                            os.path.join(data_dir, f, csv), header=1
                        )
                    }
                )

        # If file, expect to represent multiple datasets
        else:
            csv = os.path.join(data_dir, f)
            df_dict.update(read_grouped_csv(csv))

    # Plot the deconvolution
    for key, df in df_dict.items():

        interest_plot = os.path.join(output_dir, "interest_{}.png".format(key))

        plot = os.path.join(output_dir, "{}.png".format(key))

        tight_plot = os.path.join(output_dir, "tight_{}.png".format(key))

        if not os.path.exists(tight_plot):
            df.plot(
                x="X(Daltons)",
                y="Y(Counts)",
                kind="line",
                xlim=(24500, 26000),
                legend=False,
            )
            plt.ylabel("Y(Counts)")
            plt.savefig(tight_plot)
            plt.close()

        if not os.path.exists(interest_plot):
            df.plot(
                x="X(Daltons)",
                y="Y(Counts)",
                kind="line",
                xlim=(22000, 26000),
                legend=False,
            )
            plt.ylabel("Y(Counts)")
            plt.savefig(interest_plot)
            plt.close()

        if not os.path.exists(plot):
            df.plot(x="X(Daltons)", y="Y(Counts)", kind="line", legend=False)
            plt.ylabel("Y(Counts)")
            plt.savefig(plot)
            plt.close()

    # Remove blank datasets
    df_dict = remove_dataset_by_filename_content(df_dict, key_string="blank")

    # Translate keys (filenames) into ratios
    # Split across two dicitionaries to parse
    # intended ratio and expected ratio
    intended_ratio_df = ratios_from_csv(
        csv=os.path.join(args.path, "experiment_summary.csv"), df_dict=df_dict
    )

    # Remove un-needed dataset
    df_dict = remove_dataset_by_filename_content(
        df_dict, key_string="NUDT7A_p026_NU0000308a"
    )

    # Manual removal of low signal datasets
    df_dict = remove_dataset_by_filename_content(
        df_dict, key_string="Post_diffraction_crystal_mass_spectra_apr_14_16_2019_H10c"
    )
    df_dict = remove_dataset_by_filename_content(
        df_dict, key_string="Post_diffraction_crystal_mass_spectra_apr_14_16_2019_A11a"
    )
    df_dict = remove_dataset_by_filename_content(
        df_dict, key_string="Post_diffraction_crystal_mass_spectra_apr_14_16_2019_A11c"
    )
    df_dict = remove_dataset_by_filename_content(
        df_dict, key_string="Post_diffraction_crystal_mass_spectra_apr_14_16_2019_B11a"
    )
    df_dict = remove_dataset_by_filename_content(
        df_dict, key_string="Post_diffraction_crystal_mass_spectra_apr_14_16_2019_D10c"
    )
    df_dict = remove_dataset_by_filename_content(
        df_dict, key_string="Post_diffraction_crystal_mass_spectra_apr_14_16_2019_E10c"
    )
    df_dict = remove_dataset_by_filename_content(
        df_dict, key_string="Post_diffraction_crystal_mass_spectra_apr_14_16_2019_F10a"
    )
    df_dict = remove_dataset_by_filename_content(
        df_dict, key_string="Post_diffraction_crystal_mass_spectra_apr_14_16_2019_F10c"
    )
    df_dict = remove_dataset_by_filename_content(
        df_dict, key_string="Post_diffraction_crystal_mass_spectra_apr_14_16_2019_G10a"
    )
    df_dict = remove_dataset_by_filename_content(
        df_dict, key_string="Post_diffraction_crystal_mass_spectra_apr_14_16_2019_H10c"
    )

    # remove datasets that do not pass signal threshold
    # currently 1000 counts in interest area 22000-26000
    # TODO Move parameters up to non hard coded region
    df_dict = remove_datasets_with_low_signal_to_noise(df_dict)

    # From Exploratory analysis these pairs of peaks are present
    # showing that degradation into three sets of paired peaks
    expected_unlabelled_peaks = np.array([25125, 24218, 23525, 22786])
    expected_labelled_peaks = expected_unlabelled_peaks + 354

    # Process all deconvolutions to ratios of peaks
    ratio_df_list = []
    for key, df in df_dict.items():
        ratio_df = get_ratios_of_expected_peaks(
            df, expected_unlabelled_peaks, expected_labelled_peaks
        )

        # Get weights of contributing peaks based on peak height
        height_weight = (ratio_df["peak_height"] + ratio_df["labelled_peak_height"]) / (
            ratio_df["labelled_peak_height"].sum() + ratio_df["peak_height"].sum()
        )

        # Get weights of contributing peaks based on peak area
        area_weight = (ratio_df["peak_area"] + ratio_df["labelled_peak_area"]) / (
            ratio_df["labelled_peak_area"].sum() + ratio_df["peak_area"].sum()
        )

        # Get a single weighted ratio using peak heights
        weighted_height_ratio = height_weight * ratio_df["height_ratio"]
        weighted_height_ratio = weighted_height_ratio.sum()

        # Get a single weighted ratio using peak areas
        weighted_area_ratio = area_weight * ratio_df["area_ratio"]
        weighted_area_ratio = weighted_area_ratio.sum()

        # Store weights and ratios, and key in ratio dataframe
        ratio_df["height_weights"] = height_weight
        ratio_df["area_weights"] = area_weight
        ratio_df["weighted_height_ratio"] = weighted_height_ratio
        ratio_df["weighted_area_ratio"] = weighted_area_ratio
        ratio_df["intended_ratio"] = intended_ratio_df.loc[
            intended_ratio_df["f_name"] == key, "intended_ratio"
        ].iloc[0]
        ratio_df["key"] = key

        # For concatenating results
        ratio_df_list.append(ratio_df)

    ratio_df = pd.concat(ratio_df_list)

    # Split into pre crystallisation and post crystallisation
    ratio_df["pre_crystal"] = ratio_df["key"].apply(
        string_contains, match=["L_", "30min", "190506"]
    )
    pre_crystal_df = ratio_df[ratio_df["pre_crystal"] == True]
    post_crystal_df = ratio_df[ratio_df["pre_crystal"] == False]

    intended_to_weighted = dict(
        zip(pre_crystal_df["intended_ratio"], pre_crystal_df["weighted_height_ratio"])
    )

    pre_crystal_plot(pre_crystal_df)
    post_crystal_plot(post_crystal_df)

    # Add calibration from pre crystal to post crystal
    df = pre_crystal_df[["intended_ratio", "weighted_height_ratio"]]
    df = df.drop_duplicates()
    df = df.sort_values(by=["intended_ratio"])
    intended_ratio = df["intended_ratio"]
    expected_ratio = df["weighted_height_ratio"]
    interp_val = np.interp(0.75, intended_ratio, expected_ratio)

    df1 = pd.DataFrame(
        {"intended_ratio": [0.75, 0], "weighted_height_ratio": [interp_val, 0]}
    )
    df = df.append(df1)
    df = df.sort_values(by=["intended_ratio"])
    df = df.rename(columns={"weighted_height_ratio": "calibrated_ratio"})
    post_crystal_df = pd.merge(post_crystal_df, df, on="intended_ratio")

    # Split into two dataframes
    post_crystal_df["Post_diffraction"] = post_crystal_df["key"].apply(
        string_contains, match=["Post_diffraction"]
    )
    post_diffract_df = post_crystal_df[post_crystal_df["Post_diffraction"] == True]
    post_crystal_df = post_crystal_df[post_crystal_df["Post_diffraction"] == False]

    post_crystal_calibrated_plot(post_crystal_df, post_diffract_df)
    post_diffract_df.to_csv("post_diffract_data.csv")
    post_crystal_df.to_csv("post_crystal_data.csv")
