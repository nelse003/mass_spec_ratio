import os
import numpy as np
import pandas as pd
import statsmodels.api as sm
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import seaborn as sns


def plot_ratio_occupancy(
    occ_df, xlabel, f_name, occ_column_name, b_col_name="B_mean", min_cond=0.4
):
    """
    Plot occupancy vs ratio, with marker scaled by b factor

    Plot occupancy vs ratio from mass spectrscopy.
    Currently uses expected ratio.
    Fits a linear model to points with occupancy > min_condition using statsmodels.
    Plots the b factor using the size of markers.
    This is independent of the maximal and minimal B factors in the data.
    Writes plot to f_name.

    Parameters
    ----------
    occ_df: pandas.DataFrame
        Dataframe containing the occupancy,
        ratio detected in mass spectroscopy,
        and B factor

    xlabel: str
        label for the xaxis

    f_name: str
        file path for output png

    occ_column_name: str
        name of column in occ_df where occupancy is stored

    b_col_name: str
        name of the column in occ_df where b facotr is stored

    min_cond: float
        condition on occupancy where to cacluate the fit from

    Returns
    -------
    None

    Notes
    ------
    Uses

    https://stackoverflow.com/questions/47074423/how-to-get-default-blue-colour-of-matplotlib-pyplot-scatter
    https://matplotlib.org/users/dflt_style_changes.html#colors-color-cycles-and-color-maps
    """
    # Data to be plotted
    x = occ_df[occ_column_name]
    y = occ_df["Expected Ratio"]

    # Independent scale for the B factor
    # marker 10 = B factor 50
    # marker 20 = B factor 100
    # Squaring is used so the scaling is in area
    marker_area = (occ_df[b_col_name] / 5) ** 2

    # Create figure to plot onto
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)

    # Change to left and bottom axis, with zero point not at bottom of graph
    # This is for showing markers near zero
    ax.spines["bottom"].set_position("zero")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.xaxis.set_label_coords(0.5, -0.05)
    plt.xlim(0, 1.2)
    plt.ylim(-0.1, 1.1)

    plt.xlabel(xlabel)
    plt.ylabel("Ratio of labelled species")

    # Plot x=y
    xline = np.linspace(0, 1, 100)
    yline = xline
    xy_line, = ax.plot(xline, yline, "k:")

    # Set up condition for occupancy given the min_cond
    cond = (y <= 1.0) & (y >= min_cond)
    xFit = x[cond]
    yFit = y[cond]

    # # Split the markers into those that have been fitted and those that haven't
    # marker_area_fit = marker_area[cond]
    # marker_area_other = marker_area[~cond]

    # Linear fit (mx+c) to occupancy vs ratio when occ >= min_cond
    fit = sm.OLS(yFit, sm.add_constant(xFit)).fit()

    print(x)
    print(y)
    print(fit.summary())
    print(fit.params[1])
    print(np.linspace(0, 1, 100) * fit.params[1])
    print(np.linspace(0, 1, 100) * fit.params[1] + fit.params[0])

    # Plot of linear fit
    fit_line = plt.plot(
        np.linspace(0, 1, 100), np.linspace(0, 1, 100) * fit.params[1] + fit.params[0]
    )

    # Scatter plots showing the occupancy vs ratio data
    scatter_1 = ax.scatter(xFit, yFit, s=marker_area)
    # scatter_2 = ax.scatter(x[~cond], y[~cond], marker_area_other, color="C0", alpha=0.3)

    # Artist objects for use with the legend
    blue_circ = mlines.Line2D(
        [], [], color="C0", marker="o", linestyle="None", markersize=10
    )

    blue_circ_2 = mlines.Line2D(
        [], [], color="C0", marker="o", linestyle="None", markersize=20
    )

    trans_blue_circ = mlines.Line2D(
        [], [], color="C0", marker="o", linestyle="None", markersize=10, alpha=0.3
    )

    # Shrink current axis's height by 20% on the bottom
    box = ax.get_position()

    # legend usind defined artist objects and x=y line and fit_line
    # legend = plt.legend(
    #     (xy_line, blue_circ, blue_circ_2),
    #     (
    #         "Refined Occupancy\n= Ratio of Labelled Species",
    #         "B factor = 50",
    #         "B factor = 100",
    #     ),
    #     prop={"size": 8},
    #     loc="upper left",
    #     frameon=False,
    # )

    # Single point in legend
    # legend = plt.legend(
    #     (xy_line, blue_circ),
    #     (
    #         "Refined Occupancy\n= Ratio of Labelled Species",
    #         "B factor = 50",
    #     ),
    #     prop={"size": 8},
    #     loc="upper left",
    #     frameon=False,
    # )


    plt.savefig(f_name, dpi=600)
    plt.close()


def plot_occ_df(occ_df, occ_column_name, b_col_name, f_name, xlabel, metric="RSZO/OCC"):

    x = occ_df[occ_column_name]
    y = occ_df["Expected Ratio"]

    # Colour bar axis
    c = occ_df[metric]

    # Independent scale for the B factor
    # marker 10 = B factor 50
    # marker 20 = B factor 100
    # Squaring is used so the scaling is in area
    marker_area = (occ_df[b_col_name] / 5) ** 2

    # Define figure
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    ax.spines["bottom"].set_position("zero")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.xaxis.set_label_coords(0.5, -0.05)
    plt.xlim(0, 1.1)
    plt.ylim(-0.1, 1.1)

    # Plot x=y
    xline = np.linspace(0, 1, 100)
    yline = xline
    xy_line, = ax.plot(xline, yline, "k:")

    # Plot scatter
    scatter = ax.scatter(x=x, y=y, c=c, s=marker_area)

    if metric == "RSZO/OCC":
        scatter.set_clim(0, 3)
    elif metric == "RSCC":
        scatter.set_clim(0.5, 1)
    else:
        # Min / max of metric will be used
        pass

    cb = fig.colorbar(scatter, label=metric, shrink=0.75)

    blue_circ = mlines.Line2D(
        [], [], color="C0", marker="o", linestyle="None", markersize=10
    )

    blue_circ_2 = mlines.Line2D(
        [], [], color="C0", marker="o", linestyle="None", markersize=20
    )

    # Shrink current axis's height by 10% on the bottom
    box = ax.get_position()
    ax.set_position([box.x0, box.y0 + box.height * 0.1, box.width, box.height * 0.8])

    legend = ax.legend(
        (xy_line, blue_circ, blue_circ_2),
        (
            "Refined Occupancy = Ratio of Labelled Species",
            "B factor = 50",
            "B factor = 100",
        ),
        loc="upper center",
        bbox_to_anchor=(0.5, -0.15),
        prop={"size": 12},
    )

    plt.xlabel(xlabel)
    plt.ylabel("Ratio of labelled species")

    plt.savefig(f_name, dpi=600)
    plt.close()


def plot_occ_colour(occ_df, method, colour="r"):

    occ_with_std = occ_df[occ_df.mean_weighted_area_ratio.notnull()]["Occupancy"]
    ratio_with_std = occ_df[occ_df.mean_weighted_area_ratio.notnull()][
        "mean_weighted_area_ratio"
    ]
    b_with_std = occ_df[occ_df.mean_weighted_area_ratio.notnull()][
        "Average B-factor (Residue)"
    ]
    std = occ_df[occ_df.mean_weighted_area_ratio.notnull()]["std_weighted_area_ratio"]

    occ = occ_df[occ_df.mean_weighted_area_ratio.isnull()]["Occupancy"]
    ratio = occ_df[occ_df.mean_weighted_area_ratio.isnull()]["Expected Ratio"]

    # Define figure
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    ax.spines["bottom"].set_position("zero")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.xaxis.set_label_coords(0.5, -0.05)
    plt.xlim(0, 1.1)
    plt.ylim(-0.1, 1.1)

    # Plot x=y
    xline = np.linspace(0, 1, 100)
    yline = xline
    xy_line, = ax.plot(xline, yline, "k:")

    plt.errorbar(
        x=occ_with_std,
        y=ratio_with_std,
        yerr=std,
        fmt="none",
        ecolor="lightgray",
        alpha=0.5,
    )

    plt.scatter(x=occ_with_std, y=ratio_with_std, c=b_with_std)

    plt.scatter(x=occ, y=ratio, marker="d")

    plt.savefig("test_{}".format(method), dpi=600)
    plt.close()


def plot_all_regression_plots(occ_df, method_colours):

    # Define figure
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    ax.spines["bottom"].set_position("zero")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.xaxis.set_label_coords(0.5, -0.05)
    ax.set(xlim=(0, 1), ylim=(0, 1))

    xline = np.linspace(0, 1, 100)
    yline = xline
    xy_line, = ax.plot(
        xline, yline, "k:", label="Crystal occupancy\n= Ratio of labelling"
    )

    for method, method_df in occ_df.groupby("method"):

        print("AAA: {}".format(method))

        x = method_df["Occupancy"]
        y = method_df["Ratio"]

        if method == "exhaustive_search":
            print(x)
            print(y)

        if method != "refmac_superposed":
            continue

        # Remove label to plot simpler legend
        ax = sns.regplot(
            x=x,
            y=y,
            scatter=False,
            label=method.replace("_", " ").capitalize(),
            line_kws={"alpha": 0.8},
            truncate=True,
            color=method_colours[method],
        )

    plt.ylabel("Ratio of covalently labelled protein")
    plt.legend(
        (xy_line,),
        ("Crystal occupancy =\nRatio of labelling",),
        loc="upper left",
        frameon=False,
    )
    plt.tight_layout()
    plt.savefig("regression_plot_buster", dpi=600)
    plt.close()


def violin_plot_b_factor(occ_df, method_colours):

    fig = plt.figure(figsize=(8, 4.5))

    ax = sns.violinplot(
        x="method", y="Average B-factor (Residue)", data=occ_df, palette=method_colours
    )

    ax.yaxis.tick_right()
    ax.spines["top"].set_visible(False)
    ax.spines["left"].set_visible(False)
    ax.yaxis.set_label_position("right")

    method_labels = [
        method.capitalize().replace("_", "\n") for method in occ_df.method.unique()
    ]

    ax.set_xticklabels(labels=method_labels)
    plt.setp(ax.collections, alpha=0.8)
    plt.xlabel("")
    plt.ylabel("Mean B factor of covalent ligand")
    plt.savefig("b_factor_violin", dpi=600)
    plt.close()


def get_plate(key_string):

    """
    parse key to get plate in string
    """
    plate = "CI0" + key_string.split("CI0")[1].split(".d")[0].split("_")[0]
    return plate


def get_well(key_string):
    """
    Parse key yot get well in string

    Parameters
    ----------
    key_string

    Returns
    -------

    """
    well = key_string.split("CI0")[1].split(".d")[0].split("_")[1]
    if len(well) == 3:
        well = well[0] + "0" + well[1:]
    return well


def plot_delta_occ(df, df1, df2, df3, df4, df5):

    joint_df = df.merge(df1, on="crystal")

    joint_df1 = df2.merge(df3, on="crystal")

    joint_df2 = df4.merge(df5, on="crystal")

    print(joint_df.columns.values)
    ax = plt.subplot(111)
    occ_diff = joint_df["Occupancy_x"] - joint_df["Occupancy_y"]
    occ_diff_1 = joint_df1["Occupancy_x"] - joint_df1["Occupancy_y"]
    occ_diff_2 = joint_df2["Occupancy_x"] - joint_df2["Occupancy_y"]

    sns.distplot(occ_diff, rug=False, hist=False, label="phenix")
    sns.distplot(occ_diff_1, rug=False, hist=False, label="buster")
    sns.distplot(occ_diff_2, rug=False, hist=False, label="refmac")
    plt.legend()
    ax.set_xlim([0, 0.6])

    print("Phenix Mean: {} std_dev: {}".format(np.mean(occ_diff), np.std(occ_diff)))
    print("buster Mean: {} std_dev: {}".format(np.mean(occ_diff_1), np.std(occ_diff_1)))
    print("refmac mean: {} std_dev: {}".format(np.mean(occ_diff_2), np.std(occ_diff_2)))

    ocd = occ_diff.append(occ_diff_1, ignore_index=True)
    ocd = ocd.append(occ_diff_2, ignore_index=True)

    print("mean {} std {}".format(np.mean(ocd), np.std(ocd)))

    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)

    plt.xlabel("Non superposed occupancy - superposed occupancy")
    plt.ylabel("Frequency")
    plt.savefig("delta_occupancy", dpi=600)


if __name__ == "__main__":

    residue_csv = "/dls/science/groups/i04-1/elliot-dev/Work/NUDT7A_mass_spec_refinements/residue_scores.csv"
    residue_df = pd.read_csv(residue_csv)

    # exhaustive_csv = (
        "/dls/science/groups/i04-1/elliot-dev/Work/"
        "NUDT7A_mass_spec_refinements/copy_atoms/exhaustive/2019-05-29/exhaustive_minima.csv"
    )

    # b_fix
    # exhaustive_csv = (
    #     "/dls/science/groups/i04-1/elliot-dev/Work/"
    #     "NUDT7A_mass_spec_refinements/copy_atoms/exhaustive_b_fix/2019-06-17/exhaustive_minima.csv"
    # )

    exh_df = pd.read_csv(exhaustive_csv)

    exh_df = exh_df.rename(
        index=str,
        columns={"b_factor": "Average B-factor (Residue)", "occupancy": "Occupancy"},
    )
    exh_df["method"] = "exhaustive_search"

    residue_df = pd.merge(
        residue_df,
        exh_df,
        on=["crystal", "method", "Average B-factor (Residue)", "Occupancy"],
        how="outer",
    )

    mounted_df = pd.read_csv("mounted_ratios.csv")

    occ_df = pd.merge(
        residue_df, mounted_df, right_on="  Mounted Crystal ID ", left_on="crystal"
    )

    occ_df.to_csv("mounted_residue.csv")

    ratio_df = pd.read_csv("post_crystal_data.csv")
    cio_df = ratio_df[ratio_df["key"].str.contains("CI0")]
    cio_df["plate"] = cio_df.key.apply(get_plate)
    #cio_df["well"] = cio_df.key.apply(get_well)

    cio_df.to_csv("post_crystal_data_plate_well.csv")

    unique_cio_df = cio_df.drop_duplicates(subset="key")

    df_list = []
    for group, df in unique_cio_df.groupby(["plate", "intended_ratio"]):
        print(group)
        df["mean_weighted_area_ratio"] = df["weighted_area_ratio"].mean()
        df["std_weighted_area_ratio"] = df["weighted_area_ratio"].std()
        df_list.append(df)

    df_with_summary = pd.concat(df_list)
    df_with_summary_short = df_with_summary[
        [
            "plate",
            "intended_ratio",
            "mean_weighted_area_ratio",
            "std_weighted_area_ratio",
        ]
    ]
    df_with_summary_short = df_with_summary_short.drop_duplicates()
    df_with_summary_short = df_with_summary_short.rename(
        columns={"intended_ratio": "Expected Ratio"}
    )

    print(df_with_summary_short)
    print(occ_df.columns.values)

    print(occ_df["  Crystal to be Mounted "])

    occ_df["plate"] = occ_df["  Crystal to be Mounted "].apply(
        lambda s: s.split("-")[0]
    )

    occ_df = pd.merge(
        occ_df, df_with_summary_short, on=["plate", "Expected Ratio"], how="outer"
    )
    occ_df["Ratio"] = occ_df["mean_weighted_area_ratio"].fillna(
        occ_df["Expected Ratio"]
    )

    occ_df.to_csv("mounted_residue_with_error.csv")

    print(occ_df.shape)

    colours = sns.husl_palette(7)
    method_colours = {
        "phenix": colours[0],
        "phenix_superposed": colours[1],
        "buster": colours[2],
        "buster_superposed": colours[3],
        "refmac": colours[4],
        "refmac_superposed": colours[5],
        "exhaustive_search": colours[6],
    }

    violin_plot_b_factor(occ_df, method_colours)
    plot_all_regression_plots(occ_df, method_colours)

    for method, method_df in occ_df.groupby("method"):

        plot_ratio_occupancy(
            occ_df=method_df,
            f_name="LEGEND_method_occupancy_{}.png".format(method),
            xlabel="Crystallographic Occupancy".format(method),
            occ_column_name="Occupancy",
            b_col_name="Average B-factor (Residue)",
            min_cond=0,
        )
        # plot_occ_colour(occ_df=method_df, method=method)

        # plot_occ_df(
        #     occ_df=method_df,
        #     occ_column_name ="Occupancy",
        #     b_col_name="Average B-factor (Residue)",
        #     xlabel="Crystallographic Occupancy ({})".format(method),
        #     f_name="RSZO_Occ_scatter_{}.png".format(method),
        #     metric="RSZO/OCC",
        # )
        #
        # plot_occ_df(
        #     occ_df=method_df,
        #     occ_column_name ="Occupancy",
        #     b_col_name="Average B-factor (Residue)",
        #     xlabel="Crystallographic Occupancy ({})".format(method),
        #     f_name="RSCC_{}.png".format(method),
        #     metric="RSCC",
        # )

    # Difference between superposed and non superposed
    phenix_df = occ_df[occ_df["method"] == "phenix"]
    phenix_superposed_df = occ_df[occ_df["method"] == "phenix_superposed"]

    buster_df = occ_df[occ_df["method"] == "buster"]
    buster_superposed_df = occ_df[occ_df["method"] == "buster_superposed"]

    refmac_df = occ_df[occ_df["method"] == "refmac"]
    refmac_superposed_df = occ_df[occ_df["method"] == "refmac_superposed"]

    plot_delta_occ(
        phenix_df,
        phenix_superposed_df,
        buster_df,
        buster_superposed_df,
        refmac_df,
        refmac_superposed_df,
    )
