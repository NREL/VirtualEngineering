from vebio.RunFunctions import VE_params
import subprocess
import os
import numpy as np

root_path = os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir, os.pardir))
br_case_folder = os.path.dirname(__file__)
eh_case_folder = os.path.join(root_path, 'EH_OpenFOAM', 'tests', 'RushtonReact', )
eh_job_history = os.path.join(eh_case_folder, 'job_history.csv')
with open(eh_job_history) as f:
    for line in f:
        pass
    job_id = line.strip()

ve_params = VE_params.load_from_file(os.path.join(root_path, f've_params.{job_id}'), verbose=True)

rho_g = ve_params.eh_out['rho_g']
rho_x = ve_params.eh_out['rho_x']
rho_f = ve_params.eh_out['rho_f']
our_base = 100.0 # mol/m^3/hr
r = 0.75  # what is this parameter?
our_max = our_base * np.exp(-rho_f/100.0) * (r + (1.0-r) * rho_g/(rho_g+rho_x))
ko=0.01

command = f'pvpython pv_extract_our.py {our_max} {ko}'
subprocess.call(command.split())

with open(os.path.join(br_case_folder, 'our_avg.dat'), 'r') as f:
    data = [line.split('\t') for line in f.read().splitlines()]
    our_avg = float(data[-1][1])
ve_params.br_out = {'OUR': our_avg}

ve_params.write_to_file(os.path.join(root_path, f've_params.{job_id}'), verbose=True)