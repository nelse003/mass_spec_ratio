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
    plt.close()


def plot_m_z(csv_path, file_name, min_x=None, max_x=None, min_y=None, max_y=None):

    df = pd.read_csv(csv_path, skiprows=[0])

    plt.rcParams["xtick.labelsize"] = 22
    plt.rcParams["ytick.labelsize"] = 22

    df.plot(
        kind="line", x="X(Thompsons)", y="Y(Counts)", legend=False, figsize=(20, 10)
    )

    ax = plt.subplot(111)

    ax.set_xlim(left=min_x, right=max_x)
    ax.set_ylim(bottom=min_y, top=max_y)

    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)

    plt.xlabel("Mass to charge (Da/e)", fontsize=22)
    plt.ylabel("Counts", fontsize=22)
    plt.tight_layout()
    plt.savefig(dpi=600, fname=file_name)
    plt.close()


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
    plt.close()


def plot_peg_tio_zip():

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

                if folder == "/home/nelse003/peg_test/M_Z" and "0_05" in f:
                    max_y = 0.08
                    min_y = 0.00
                else:
                    max_y = None
                    min_y = None

                plot_m_z(
                    csv_path=os.path.join(folder, f),
                    file_name=os.path.join(
                        folder, "m_z_{}".format(f.replace(".CSV", ".png"))
                    ),
                    min_x=600,
                    max_y=max_y,
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


if __name__ == "__main__":

    chromatogram_folder = "/home/enelson/uHPLC_method/chromatogram"

    for f in os.listdir(chromatogram_folder):
        if ".CSV" in f:
            plot_chromatogram(
                csv_path=os.path.join(chromatogram_folder, f),
                file_name=os.path.join(
                    chromatogram_folder,
                    "Chromatogram_{}".format(f.replace(".CSV", ".png")),
                ),
            )

    m_z_folder = "/home/enelson/uHPLC_method/m_z"

    for f in os.listdir(m_z_folder):
        if ".CSV" in f:
            plot_m_z(
                csv_path=os.path.join(m_z_folder, f),
                file_name=os.path.join(
                    m_z_folder, "m_z_{}".format(f.replace(".CSV", ".png"))
                ),
                min_x=600,
                max_y=None,
            )

    deconvolution_folder = "/home/enelson/uHPLC_method/deconvolution"

    for f in os.listdir(deconvolution_folder):
        if ".CSV" in f:
            plot_deconvolution(
                csv_path=os.path.join(deconvolution_folder, f),
                file_name=os.path.join(
                    deconvolution_folder, "deconv_{}".format(f.replace(".CSV", ".png"))
                ),
            )

            plot_deconvolution(
                csv_path=os.path.join(deconvolution_folder, f),
                file_name=os.path.join(
                    deconvolution_folder,
                    "tight_deconv_{}".format(f.replace(".CSV", ".png")),
                ),
                min_x=25000,
                max_x=26000,
            )

    plot_deconvolution(
        csv_path="/home/enelson/ion_suppression/NUDT7A_p026_cps086_12.1.CSV",
        file_name="/home/enelson/ion_suppression/deconv_NUDT7A_p026_cps086_12.png",
    )

    plot_deconvolution(
        csv_path="/home/enelson/ion_suppression/NUDT7A_p026_cps086_12.1.CSV",
        file_name="/home/enelson/ion_suppression/tight_deconv_NUDT7A_p026_cps086_12.png",
        min_x=22000,
        max_x=26000,
    )
