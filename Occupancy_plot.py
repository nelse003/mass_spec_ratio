import os
import numpy as np
import pandas as pd
import statsmodels.api as sm
import matplotlib.pyplot as plt
import matplotlib.lines as mlines

plt.rcParams.update({"font.size": 18})
plt.rcParams["figure.figsize"] = [10, 10]


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
    plt.xlim(0, 1.0)
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

    # Split the markers into those that have been fitted and those that haven't
    marker_area_fit = marker_area[cond]
    marker_area_other = marker_area[~cond]

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
    scatter_1 = ax.scatter(xFit, yFit, s=marker_area_fit)
    scatter_2 = ax.scatter(x[~cond], y[~cond], marker_area_other, color="C0", alpha=0.3)

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

    # legend usind defined artist objects and x=y line and fit_line
    legend = plt.legend(
        (xy_line, trans_blue_circ, fit_line, blue_circ, blue_circ_2),
        (
            "Refined Occupancy = Ratio of Labelled Species",
            "Data not fitted",
            "Fit to ratio >= {}\n {:.2f}x{:+.2f}. \n R-squared {:.3f}".format(
                min_cond, fit.params[1], fit.params[0], fit.rsquared
            ),
            "B factor = 50",
            "B factor = 100",
        ),
        prop={"size": 12},
    )

    plt.savefig(f_name, dpi=600)


if __name__ == "__main__":

    ref_dir = "/dls/science/groups/i04-1/elliot-dev/Work/NUDT7A_mass_spec_refinements/copy_atoms_190525_buster"

    occupancies = {}
    for refinement_folder in os.listdir(ref_dir):
        if os.path.isdir(os.path.join(ref_dir, refinement_folder)):

            bound_state = os.path.join(
                ref_dir, refinement_folder, "refine.split.bound-state.pdb"
            )

            occ_file =  os.path.join(ref_dir, refinement_folder, "occ.txt")

            if os.path.isfile(bound_state):
                os.system("ccp4-python ccp4/occ_b.py {} {} E 1 ".format(bound_state,occ_file))
                with open(occ_file) as f:
                    mean_occ, mean_b, std_b =  f.readline().split(',')
                    occupancies[refinement_folder] = (float(mean_occ), float(mean_b), float(std_b))


    occ_df = pd.DataFrame.from_dict(occupancies, orient='index',columns=["occupancy", "b_mean", "b_std"])
    occ_df['crystal'] = occ_df.index

    occ_correct = "occ_correct.csv"
    occ_correct_df = pd.read_csv(occ_correct)
    mounted_df = pd.read_csv("mounted_ratios.csv")

    occ_df = pd.merge(
        occ_df, mounted_df, right_on="  Mounted Crystal ID ", left_on="crystal"
    )

    plot_ratio_occupancy(
        occ_df=occ_df,
        f_name="Occupancy_B_factor_buster.png",
        xlabel="Crystallographic Occupancy (buster)",
        occ_column_name="occupancy",
        b_col_name="b_mean"
    )
