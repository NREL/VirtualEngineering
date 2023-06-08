import os
from vebio.FileModifiers import write_file_with_replacements
from vebio.RunFunctions import VE_params

root_path = os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir, os.pardir))
eh_case_folder = os.path.join(root_path, 'EH_OpenFOAM', 'tests', 'RushtonReact', )

eh_job_history = os.path.join(eh_case_folder, 'job_history.csv')
with open(eh_job_history) as f:
    for line in f:
        pass
    job_id = line

ve = VE_params.load_from_file(os.path.join(root_path, 've_params.{job_id}'), verbose=True)
# Make changes to the fvOptions file based on replacement options
fvOptions = {}
fvOptions['rho_g'] = ve.eh_out['rho_g']
fvOptions['rho_x'] = ve.eh_out['rho_x']
fvOptions['rho_f'] = ve.eh_out['rho_f']
write_file_with_replacements(os.path.join(os.path.dirname(__file__), 'constant', 'fvOptions'), fvOptions)

# Make changes to the controlDict file based on replacement options
controlDict = {}
controlDict['endTime'] = ve.br_in['t_final']
write_file_with_replacements(os.path.join(os.path.dirname(__file__), 'system', 'controlDict'), controlDict)

