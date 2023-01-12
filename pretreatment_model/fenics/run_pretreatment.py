from Pretreatment import Pretreatment
import cProfile
import time

def main():

    t1 = time.time()

    PT = Pretreatment()

    PT.generate_mesh(nn=16)

    PT.build_functions(degree=2)

    PT.build_problem()

    PT.solve(t_final=1200.0, dt=2.0, save_every=60.0)

    t2 = time.time()

    print(f'Pretreatment wallclock time: {t2-t1:.2f} s.')


if __name__ == "__main__":
    profiler = cProfile.Profile()
    profiler.enable()
    main()
    profiler.disable()
    profiler.dump_stats("profiling.prof")
