import os
import matplotlib.pyplot as plt
import pandas as pd


def plot_chromatogram(csv_path, file_name):
    """Plot chromatogram from csv"""

    df = pd.read_csv(csv_path, skiprows=[0])

    plt.rcParams["xtick.labelsize"] = 22
    plt.rcParams["ytick.labelsize"] = 22

    df.plot(kind="line", x="X(Minutes)", y="Y(Counts)", legend=False, figsize=(20, 10))

    ax = plt.subplot(111)

    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)

    plt.xlabel("Aquistion time (Minutes)", fontsize=22)
    plt.ylabel("Counts", fontsize=22)
    plt.savefig(dpi=600, fname=file_name)


def plot_m_z(csv_path, file_name, min_x=None, max_x=None):

    df = pd.read_csv(csv_path, skiprows=[0])

    plt.rcParams["xtick.labelsize"] = 22
    plt.rcParams["ytick.labelsize"] = 22

    df.plot(
        kind="line", x="X(Thompsons)", y="Y(Counts)", legend=False, figsize=(20, 10)
    )

    ax = plt.subplot(111)

    ax.set_xlim(left=min_x, right=max_x)

    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)

    plt.xlabel("Mass to charge (Da/e)", fontsize=22)
    plt.ylabel("Counts", fontsize=22)
    plt.savefig(dpi=600, fname=file_name)


def plot_deconvolution(csv_path, file_name, min_x=None, max_x=None):

    df = pd.read_csv(csv_path, skiprows=[0])

    plt.rcParams["xtick.labelsize"] = 22
    plt.rcParams["ytick.labelsize"] = 22

    df.plot(kind="line", x="X(Daltons)", y="Y(Counts)", legend=False, figsize=(20, 10))

    ax = plt.subplot(111)

    ax.set_xlim(left=min_x, right=max_x)

    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)

    plt.xlabel("Deconvoluted mass (Da)", fontsize=22)
    plt.ylabel("Counts", fontsize=22)
    plt.savefig(dpi=600, fname=file_name)


if __name__ == "__main__":

    peg_test_chromatograms = "/home/nelse003/peg_test/Chromatograms"
    tio_chromatograms = "/home/nelse003/tio_removal/chromatograms"
    zip_tip_chromatograms = "/home/nelse003/zip_tip/Chromatograms"

    chromatograms = [peg_test_chromatograms, tio_chromatograms, zip_tip_chromatograms]

    for folder in chromatograms:
        for f in os.listdir(folder):
            if ".CSV" in f:
                plot_chromatogram(
                    csv_path=os.path.join(folder, f),
                    file_name=os.path.join(
                        folder, "Chromatogram_{}".format(f.replace(".CSV", ".png"))
                    ),
                )

    peg_test_m_z = "/home/nelse003/peg_test/M_Z"
    tio_m_z = "/home/nelse003/tio_removal/m_z"
    zip_tip_m_z = "/home/nelse003/zip_tip/m_z"

    m_z = [peg_test_m_z, tio_m_z, zip_tip_m_z]

    for folder in m_z:
        for f in os.listdir(folder):
            if ".CSV" in f:
                plot_m_z(
                    csv_path=os.path.join(folder, f),
                    file_name=os.path.join(
                        folder, "m_z_{}".format(f.replace(".CSV", ".png"))
                    ),
                    min_x=600,
                )

    peg_test_deconv = "/home/nelse003/peg_test/deconvolutions"
    tio_deconv = "/home/nelse003/tio_removal/deconvolutions"
    zip_tip_deconv = "/home/nelse003/zip_tip/deconvolutions"

    deconv = [peg_test_deconv, tio_deconv, zip_tip_deconv]

    for folder in deconv:
        for f in os.listdir(folder):
            if ".CSV" in f:

                plot_deconvolution(
                    csv_path=os.path.join(folder, f),
                    file_name=os.path.join(
                        folder, "deconv_{}".format(f.replace(".CSV", ".png"))
                    ),
                )

                plot_deconvolution(
                    csv_path=os.path.join(folder, f),
                    file_name=os.path.join(
                        folder, "tight_deconv_{}".format(f.replace(".CSV", ".png"))
                    ),
                    min_x=25000,
                    max_x=26000,
                )
