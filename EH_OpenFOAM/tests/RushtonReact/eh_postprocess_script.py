# if EH ran successfully 
import os
import sys
import numpy as np
from vebio.RunFunctions import VE_params

root_path = os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir, os.pardir, os.pardir))
job_id = sys.argv[1]
ve = VE_params.load_from_file(os.path.join(root_path, f've_params.{job_id}'), verbose=True)
dilution_factor = ve.eh_in['fis_0']/ve.pt_out['fis_0']

integrated_quantities = np.genfromtxt('integrated_quantities.dat') # mol/L

ve.eh_out['rho_g'] = float(integrated_quantities[-1, -3])
ve.eh_out['rho_x'] = float(integrated_quantities[-1, -2])
ve.eh_out['rho_sl'] = float(integrated_quantities[-1, -1])
ve.eh_out['rho_f'] = float(ve.pt_out['rho_f']*dilution_factor)

ve.write_to_file(os.path.join(root_path, f've_params.{job_id}'), verbose=True)