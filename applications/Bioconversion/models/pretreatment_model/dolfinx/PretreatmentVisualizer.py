import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import os


def get_plot_color(cmap, tt_frac):

    return cmap(0.9 * tt_frac)


def draw_figures(
    t_final,
    path_to_data_files,
    tt_list=[0, 30, 300, 600, 1200],
    plot_fortran_truth=True,
):

    cmap = cm.get_cmap("viridis")

    fig1, ax1 = plt.subplots(2, 2, figsize=(8, 6), dpi=300)
    fig2, ax2 = plt.subplots(2, 2, figsize=(8, 6), dpi=300)

    ax1[1, 0].set_xlabel("Distance (cm)")
    ax1[1, 1].set_xlabel("Distance (cm)")
    ax1[0, 0].set_ylabel("Steam (M)")
    ax1[0, 1].set_ylabel("Liquid Fraction")
    ax1[1, 0].set_ylabel("Temperature (K)")
    ax1[1, 1].set_ylabel("Acid (M)")

    ax2[1, 0].set_xlabel("Distance (cm)")
    ax2[1, 1].set_xlabel("Distance (cm)")
    ax2[0, 0].set_ylabel("Xylan Fraction")
    ax2[0, 1].set_ylabel("Xylooligomers (M)")
    ax2[1, 0].set_ylabel("Xylose (M)")
    ax2[1, 1].set_ylabel("Furfural (mM)")

    for tt in tt_list:
        data_filename = os.path.join(path_to_data_files, f"data_t{tt:05.0f}s.csv")
        data = np.genfromtxt(data_filename, delimiter=",", skip_header=1)

        include_reflected = True

        if include_reflected:
            data_reflected = np.copy(data)
            data_reflected = np.flipud(data_reflected)
            data_reflected[:, 0] *= -1
            data = np.vstack((data_reflected[:-1, :], data))

        color = get_plot_color(cmap, tt / t_final)
        mew = 0.5
        marker = "."

        ax1[0, 0].plot(
            data[:, 0],
            data[:, 1],
            marker,
            markeredgewidth=mew,
            markersize=2,
            label=f"{tt:.0f} s",
            color=color,
        )
        # Liquid Fraction
        ax1[0, 1].plot(
            data[:, 0],
            data[:, 2],
            marker,
            markeredgewidth=mew,
            markersize=2,
            color=color,
        )
        # Temperature
        ax1[1, 0].plot(
            data[:, 0],
            data[:, 3],
            marker,
            markeredgewidth=mew,
            markersize=2,
            color=color,
        )
        # Acid
        ax1[1, 1].plot(
            data[:, 0],
            data[:, 4],
            marker,
            markeredgewidth=mew,
            markersize=2,
            color=color,
        )

        # Xylan Fraction
        ax2[0, 0].plot(
            data[:, 0],
            data[:, 5],
            marker,
            markeredgewidth=mew,
            markersize=2,
            label=f"{tt:.0f} s",
            color=color,
        )
        # Xylooligomer
        ax2[0, 1].plot(
            data[:, 0],
            data[:, 6],
            marker,
            markeredgewidth=mew,
            markersize=2,
            color=color,
        )
        # Xylose
        ax2[1, 0].plot(
            data[:, 0],
            data[:, 7],
            marker,
            markeredgewidth=mew,
            markersize=2,
            color=color,
        )
        # Furfural
        ax2[1, 1].plot(
            data[:, 0],
            data[:, 8],
            marker,
            markeredgewidth=mew,
            markersize=2,
            color=color,
        )

    # plot truth curves:
    if plot_fortran_truth:
        for z, truth_file in enumerate(["truth"]):
            for ctr_truth, tt_truth in enumerate([30, 300, 600, 1200]):
                try:
                    truth_filename = os.path.join(
                        "truth", f"{truth_file}_{tt_truth}s.dat"
                    )
                    truth_data = np.genfromtxt(truth_filename, skip_header=1)
                except:
                    print("Could not find truth data file to read.")
                    break

                if include_reflected:
                    truth_data_reflected = np.copy(truth_data)
                    truth_data_reflected = np.flipud(truth_data_reflected)
                    truth_data_reflected[:, 0] *= -1
                    truth_data = np.vstack((truth_data_reflected[:-1, :], truth_data))

                rho_x = 730.0
                rho_os = 1400.0
                f_x0 = 0.26
                eps_p0 = 0.8
                c_acid0 = 0.1 * 1e3
                eps_l0 = 0.25

                aa = rho_x / rho_os
                bb = f_x0 / (1.0 - f_x0) + aa
                cc = 1.0 - eps_p0

                ss1 = aa**2 * cc**2
                ss2 = 2.0 * (aa - 2.0) * bb * cc * truth_data[:, 4]
                ss3 = bb**2 * truth_data[:, 4] ** 2

                ss = np.sqrt(ss1 - ss2 + ss3)

                np_eps_p = -(aa * cc + bb * (truth_data[:, 4] - 2.0) + ss) / (2.0 * bb)

                ls_list = ["-", "--"]
                ls = ls_list[z]
                lw = 0.75

                color = get_plot_color(cmap, tt_truth / t_final)

                ax1[0, 0].plot(
                    truth_data[:, 0],
                    1e3 * truth_data[:, 1],
                    linestyle=ls,
                    linewidth=lw,
                    color=color,
                    label=f"{truth_file} {tt_truth} s",
                )
                ax1[0, 1].plot(
                    truth_data[:, 0],
                    truth_data[:, 2],
                    linestyle=ls,
                    linewidth=lw,
                    color=color,
                )
                ax1[1, 0].plot(
                    truth_data[:, 0],
                    truth_data[:, 3],
                    linestyle=ls,
                    linewidth=lw,
                    color=color,
                )
                ax1[1, 1].plot(
                    truth_data[:, 0],
                    1e-3 * c_acid0 * eps_l0 / truth_data[:, 2],
                    linestyle=ls,
                    linewidth=lw,
                    color=color,
                )

                ax2[0, 0].plot(
                    truth_data[:, 0],
                    truth_data[:, 4] / (1 - np_eps_p),
                    linestyle=ls,
                    linewidth=lw,
                    color=color,
                    label=f"{truth_file} {tt_truth} s",
                )
                ax2[0, 1].plot(
                    truth_data[:, 0],
                    1e3 * truth_data[:, 5] / truth_data[:, 2],
                    linestyle=ls,
                    linewidth=lw,
                    color=color,
                )
                ax2[1, 0].plot(
                    truth_data[:, 0],
                    1e3 * truth_data[:, 6] / truth_data[:, 2],
                    linestyle=ls,
                    linewidth=lw,
                    color=color,
                )
                ax2[1, 1].plot(
                    truth_data[:, 0],
                    1e6 * truth_data[:, 7] / truth_data[:, 2],
                    linestyle=ls,
                    linewidth=lw,
                    color=color,
                )

    fig1.legend(loc=7, bbox_to_anchor=(1.22, 0.5))
    fig1.tight_layout()
    fig1.savefig("output/figure_1.png", bbox_inches="tight")
    fig1.savefig("output/figure_1.pdf", bbox_inches="tight")

    fig2.legend(loc=7, bbox_to_anchor=(1.22, 0.5))
    fig2.tight_layout()
    fig2.savefig("output/figure_2.png", bbox_inches="tight")
    fig2.savefig("output/figure_2.pdf", bbox_inches="tight")
