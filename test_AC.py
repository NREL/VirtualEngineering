from ipywidgets import *
import os
import numpy as np

# imports from vebio modules
from vebio.WidgetFunctions import csv2widget_collection
from vebio.RunFunctions import Pretreatment, Feedstock, EnzymaticHydrolysis, Bioreactor, VE_params
from vebio.OptimizationFunctions import Optimization


hpc_run = True

# Set Virtual Engineering Options
fs_options = csv2widget_collection("feedstock_params.csv")
pt_options = csv2widget_collection("pretreatment_params.csv")
eh_options = csv2widget_collection("enzymatic_hydrolysis_params.csv")
br_options = csv2widget_collection("bioreactor_params.csv")

obj_widget = widgets.Dropdown(
    options=[('Biorector:            OUR',   ('br_out', 'our'))],value=('br_out', 'our'))

eh_options.fis_0.is_control.value = True  

Opt_lf = Optimization(fs_options, pt_options, eh_options, br_options, obj_widget, hpc_run)
Opt_lf.EH_model.model_type = 'CFD Surrogate'
Opt_lf.BR_model.model_type = 'CFD Surrogate'
print('x vector order', Opt_lf.var_names)

Opt_hf = Optimization(fs_options, pt_options, eh_options, br_options, obj_widget, hpc_run)
Opt_hf.EH_model.model_type = 'CFD Simulation'
Opt_hf.BR_model.model_type = 'CFD Simulation'
print('x vector order', Opt_hf.var_names)

Opt_mix = Optimization(fs_options, pt_options, eh_options, br_options, obj_widget, hpc_run)
Opt_mix.EH_model.model_type = 'CFD Surrogate'
Opt_mix.BR_model.model_type = 'CFD Simulation'
print('x vector order', Opt_mix.var_names)


def lf_simulation(x):
    return Opt_lf.run_models_with_new_values(x)

def hf_simulation(x):
    Opt_hf.run_models_with_new_values(x)
    root_path = os.path.join(os.path.dirname(__file__))
    eh_case_folder = os.path.join(root_path, 'EH_OpenFOAM', 'tests', 'RushtonReact', )

    eh_job_history = os.path.join(eh_case_folder, 'job_history.csv')
    with open(eh_job_history) as f:
        for line in f:
            pass
        job_id = line.strip()

    ve = VE_params.load_from_file(os.path.join(root_path, f've_params.{job_id}'), verbose=True)
    return ve.br_out['OUR']

def mix_simulation(x):
    Opt_mix.run_models_with_new_values(x)
    root_path = os.path.join(os.path.dirname(__file__))
    eh_case_folder = os.path.join(root_path, 'EH_OpenFOAM', 'tests', 'RushtonReact', )

    eh_job_history = os.path.join(eh_case_folder, 'job_history.csv')
    with open(eh_job_history) as f:
        for line in f:
            pass
        job_id = line.strip()

    ve = VE_params.load_from_file(os.path.join(root_path, f've_params.{job_id}'), verbose=True)
    return ve.br_out['OUR']


['initial_acid_conc', 'lambda_e', 'fis_0']
x = [0.0001, 30, 0.05]
our = lf_simulation(x)
print('Low Fidelity simulation OUR:', our)

# our = mix_simulation(x)
# print('Mixed Fidelity simulation OUR:', our)
