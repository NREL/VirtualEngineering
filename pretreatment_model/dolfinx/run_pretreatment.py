from Pretreatment import Pretreatment
from PretreatmentVisualizer import draw_figures
import cProfile
import time

import os
import sys
root_path = os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir, os.pardir))
sys.path.append(root_path)

from vebio.RunFunctions import VE_params


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
        print(show_plots)
        draw_figures(t_final, PT.path_to_data_files)

    t2 = time.time()

    if verbose:
        print(f"Pretreatment wallclock time: {t2-t1:.2f} s.")

    output_dict = PT.integrated_quantities.copy()

    return output_dict


if __name__ == "__main__":
    if len(sys.argv) == 1:
        print('Running Preatreatment with defualt parameters')
        profiler = cProfile.Profile()
        profiler.enable()
        run_pt()
        profiler.disable()
        profiler.dump_stats("profiling.prof")
    elif len(sys.argv) == 3:
        ve_params = VE_params.load_from_file('ve_params.yml')
        output_dict = run_pt(ve_params=ve_params, verbose=sys.argv[1], show_plots=(sys.argv[2] == 'True'))       
        ve_params.pt_out = output_dict
        ve_params.write_to_file('ve_params.yml')
    else:
        print('Need 2 arguments: Defualt are verbose=True, show_plots=True')
