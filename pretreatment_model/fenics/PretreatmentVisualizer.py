import matplotlib.pyplot as plt
from fenics import FunctionSpace, project
import numpy as np

def update_figure_1(PT, tt, line_id):
    if not hasattr(PT, "ax1"):
        PT.fig1, PT.ax1 = plt.subplots(2, 2, figsize=(8, 6), dpi=300)

        PT.ax1[0, 0].set_ylabel("Steam (M)")

        PT.ax1[0, 1].set_ylabel("Liquid Fraction")

        PT.ax1[1, 0].set_xlabel("Distance (cm)")
        PT.ax1[1, 0].set_ylabel("Temperature (K)")

        PT.ax1[1, 1].set_xlabel("Distance (cm)")
        PT.ax1[1, 1].set_ylabel("Acid (M)")

        PT.color_arr = ["k", "C0", "C1", "C2", "C3"]
        PT.linestyle_arr = ["-", "-", "-", "-", "-"]

        PT.dpS = FunctionSpace(PT.mesh, 'DP', 0)

    color = PT.color_arr[line_id]
    PT.mew = 0.5

    # Steam
    c_s_vec = project(PT.u[PT.c_s_id], PT.S, solver_type="cg")
    c_s_vec = c_s_vec.vector()[:] / 1000.0
    c_s_vec = np.hstack((c_s_vec, c_s_vec))
    xx = PT.S.tabulate_dof_coordinates().flatten()
    xpts = 100.0*np.hstack((-xx, xx))
    xpts_id = np.argsort(xpts)
    PT.ax1[0, 0].plot(xpts[xpts_id], c_s_vec[xpts_id], "x", markeredgewidth=PT.mew, markersize=2, label=f"{tt:.1f} s", color=color)

    # Liquid Fraction
    eps_l_vec = project(PT.u[PT.eps_l_id], PT.S, solver_type="cg")
    eps_l_vec = eps_l_vec.vector()[:]
    eps_l_vec = np.hstack((eps_l_vec, eps_l_vec))
    PT.ax1[0, 1].plot(
        xpts[xpts_id], eps_l_vec[xpts_id], "x", markeredgewidth=PT.mew, markersize=2, color=color
    )

    # Temperature
    T_vec = project(PT.u[PT.T_id], PT.S, solver_type="cg")
    T_vec = T_vec.vector()[:]
    T_vec = np.hstack((T_vec, T_vec))
    PT.ax1[1, 0].plot(xpts[xpts_id], T_vec[xpts_id], "x", markeredgewidth=PT.mew, markersize=2, color=color)

    # Acid
    acid_vec = project(PT.c_acid, PT.S, solver_type="cg")
    acid_vec = acid_vec.vector()[:] / 1000.0
    acid_vec = np.hstack((acid_vec, acid_vec))
    PT.ax1[1, 1].plot(xpts[xpts_id], acid_vec[xpts_id], "x", markeredgewidth=PT.mew, markersize=2, color=color)

    return PT

def update_figure_2(PT, tt, line_id):
    if not hasattr(PT, "ax2"):
        PT.fig2, PT.ax2 = plt.subplots(2, 2, figsize=(8, 6), dpi=300)

        PT.ax2[0, 0].set_ylabel("Xylan Fraction")

        PT.ax2[0, 1].set_ylabel("Xylooligomers (M)")

        PT.ax2[1, 0].set_xlabel("Distance (cm)")
        PT.ax2[1, 0].set_ylabel("Xylose (M)")

        PT.ax2[1, 1].set_xlabel("Distance (cm)")
        PT.ax2[1, 1].set_ylabel("Furfural (mM)")

    color = PT.color_arr[line_id]

    # Xylan
    f_x_vec = project(PT.u[PT.fstar_x_id] / (1.0 - PT.eps_p), PT.S, solver_type="cg")
    f_x_vec = f_x_vec.vector()[:]
    f_x_vec = np.hstack((f_x_vec, f_x_vec))
    xx = PT.S.tabulate_dof_coordinates().flatten()
    xpts = 100.0*np.hstack((-xx, xx))
    xpts_id = np.argsort(xpts)
    PT.ax2[0, 0].plot(xpts[xpts_id], f_x_vec[xpts_id], "x", markeredgewidth=PT.mew, markersize=2, label=f"{tt:.1f} s", color=color)

    # Xylooligomers
    c_xo_vec = project(PT.u[PT.cstar_xo_id] / PT.u[PT.eps_l_id], PT.S, solver_type="cg")
    c_xo_vec = c_xo_vec.vector()[:] / 1000.0
    c_xo_vec = np.hstack((c_xo_vec, c_xo_vec))
    PT.ax2[0, 1].plot(xpts[xpts_id], c_xo_vec[xpts_id], "x", markeredgewidth=PT.mew, markersize=2, color=color)

    # Xylose
    c_xy_vec = project(PT.u[PT.cstar_xy_id] / PT.u[PT.eps_l_id], PT.S, solver_type="cg")
    c_xy_vec = c_xy_vec.vector()[:] / 1000.0
    c_xy_vec = np.hstack((c_xy_vec, c_xy_vec))
    PT.ax2[1, 0].plot(xpts[xpts_id], c_xy_vec[xpts_id], "x", markeredgewidth=PT.mew, markersize=2, color=color)

    # Furfural
    c_f_vec = project(PT.u[PT.cstar_f_id] / PT.u[PT.eps_l_id], PT.S, solver_type="cg")
    c_f_vec = c_f_vec.vector()[:]
    c_f_vec = np.hstack((c_f_vec, c_f_vec))
    PT.ax2[1, 1].plot(xpts[xpts_id], c_f_vec[xpts_id], "x", markeredgewidth=PT.mew, markersize=2, color=color)

    return PT

def finalize_figures(PT):
    # x(cm)   Steam(molcm3)   Liquidfrac   Temperature(K)   Xylanfrac   Xylog(molcm3)   Xylose(molcm3)   Furfural(molcm3)


    # plot truth curves:
    for ctr_truth, tt_truth in enumerate([30, 300, 600, 1200]):

        for z, truth_file in enumerate(["truth", "old_truth"]):
            try:
                truth_data = np.genfromtxt(f"{truth_file}_{tt_truth}s.dat", skip_header=1)
            except:
                print(f'Could not find "{truth_file}_{tt_truth}s.dat" file to read.')
                break


            aa = PT.rho_x/PT.rho_os
            bb = PT.f_x0/(1.0-PT.f_x0) + aa
            cc = 1.0 - PT.eps_p0

            ss1 = aa**2 * cc**2
            ss2 = 2.0 * (aa - 2.0) * bb * cc * truth_data[:, 4]
            ss3 = bb**2 * truth_data[:, 4]**2

            ss = np.sqrt(ss1 - ss2 + ss3)

            np_eps_p = -(aa*cc + bb*(truth_data[:, 4] - 2.0) + ss)/ (2.0*bb)

            ls_list = ['-', '--']
            ls = ls_list[z]
            lw = 1.0


            PT.ax1[0, 0].plot(truth_data[:, 0], 1e3*truth_data[:, 1], linestyle=ls, linewidth=lw, color=f"C{ctr_truth}", label=f'{truth_file} {tt_truth}.0s')
            PT.ax1[0, 1].plot(truth_data[:, 0], truth_data[:, 2], linestyle=ls, linewidth=lw, color=f"C{ctr_truth}")
            PT.ax1[1, 0].plot(truth_data[:, 0], truth_data[:, 3], linestyle=ls, linewidth=lw, color=f"C{ctr_truth}")
            PT.ax1[1, 1].plot(truth_data[:, 0], 1e-3*PT.c_acid0*PT.eps_l0/truth_data[:, 2], linestyle=ls, linewidth=lw, color=f"C{ctr_truth}")

            PT.ax2[0, 0].plot(truth_data[:, 0], truth_data[:, 4]/(1-np_eps_p), linestyle=ls, linewidth=lw, color=f"C{ctr_truth}", label=f'{truth_file} {tt_truth}.0s')
            PT.ax2[0, 1].plot(truth_data[:, 0], 1e3*truth_data[:, 5]/truth_data[:, 2], linestyle=ls, linewidth=lw, color=f"C{ctr_truth}")
            PT.ax2[1, 0].plot(truth_data[:, 0], 1e3*truth_data[:, 6]/truth_data[:, 2], linestyle=ls, linewidth=lw, color=f"C{ctr_truth}")
            PT.ax2[1, 1].plot(truth_data[:, 0], 1e6*truth_data[:, 7]/truth_data[:, 2], linestyle=ls, linewidth=lw, color=f"C{ctr_truth}")

    # PT.ax1[0, 0].set_ylim(0, 0.16)
    # PT.ax1[0, 1].set_ylim(0.24, 0.41)
    # PT.ax1[1, 0].set_ylim(360, 470)
    # PT.ax1[1, 1].set_ylim(0.03, 0.18)
    PT.fig1.legend(loc=7, bbox_to_anchor=(1.22, 0.5))


    PT.fig1.tight_layout()
    PT.fig1.savefig("output/figure_1.png", bbox_inches="tight")
    PT.fig1.savefig("output/figure_1_smooth.pdf", bbox_inches="tight")


    # PT.ax2[0, 0].set_ylim(0.21, .27)
    # PT.ax2[0, 1].set_ylim(0, 0.025)
    # PT.ax2[1, 0].set_ylim(0, 0.2)
    # PT.ax2[1, 1].set_ylim(0, 0.4)
    PT.fig2.legend(loc=7, bbox_to_anchor=(1.22, 0.5))

    PT.fig2.tight_layout()
    PT.fig2.savefig("output/figure_2.png", bbox_inches="tight")
    PT.fig2.savefig("output/figure_2_smooth.pdf", bbox_inches="tight")

    u_list = PT.u.split(deepcopy=True)

    temp = []

    for val in u_list:
        temp.append(val.vector()[:])

    temp = np.array(temp)

    # np.savetxt('truth.csv', temp, delimiter=',')

    # truth = np.genfromtxt('truth.csv', delimiter=',')

    # abs_err = np.abs(temp - truth)

    # if np.any(abs_err > 1e-9):
    #     print('ABS ERROR.')
    # else:
    #     print("Test passed.")

    # print(f'Max error: {np.amax(abs_err)}')

    return PT


