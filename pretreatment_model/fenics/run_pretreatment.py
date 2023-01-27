from Pretreatment import Pretreatment
import cProfile
import time

def run_pt(ve_params=None, verbose=True, show_plots=True):

    t1 = time.time()

    PT = Pretreatment(verbose, show_plots)

    PT.generate_mesh(nn=16)

    PT.build_functions(degree=2)

    PT.build_problem(ve_params)

    t_final = ve_params.pt_in['final_time'] # in seconds
    if ve_params is None:
        t_final = 1200.0   
    PT.solve(t_final=t_final, dt=2.0, save_every=60.0)

    t2 = time.time()

    print(f'Pretreatment wallclock time: {t2-t1:.2f} s.')

    output_dict = {}
    return output_dict


if __name__ == "__main__":
    profiler = cProfile.Profile()
    profiler.enable()
    run_pt()
    profiler.disable()
    profiler.dump_stats("profiling.prof")
