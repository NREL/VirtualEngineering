from Pretreatment import Pretreatment
from PretreatmentVisualizer import draw_figures
import cProfile
import time


def run_pt(ve_params=None, verbose=True, show_plots=True):

    t1 = time.time()

    PT = Pretreatment(verbose, show_plots)

    PT.generate_mesh(nn=32)

    PT.build_functions(degree=2)

    PT.build_problem(ve_params)

    if ve_params is None:
        t_final = 1200.0
    else:
        t_final = ve_params.pt_in["final_time"]  # in seconds

    PT.solve(t_final=t_final, dt=2.0, save_every=30.0)

    if show_plots:
        draw_figures(t_final, PT.path_to_data_files)

    t2 = time.time()

    if verbose:
        print(f"Pretreatment wallclock time: {t2-t1:.2f} s.")

    output_dict = PT.integrated_quantities.copy()

    return output_dict


if __name__ == "__main__":
    profiler = cProfile.Profile()
    profiler.enable()
    run_pt()
    profiler.disable()
    profiler.dump_stats("profiling.prof")
