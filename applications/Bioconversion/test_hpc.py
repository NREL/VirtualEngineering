from ipywidgets import *
import os
import numpy as np

# imports from vebio modules
from vebio.WidgetFunctions import csv2widget_collection
from vebio.RunFunctions import Pretreatment, Feedstock, EnzymaticHydrolysis, Bioreactor
from vebio.OptimizationFunctions import Optimization


hpc_run = True

# Set Virtual Engineering Options
fs_options = csv2widget_collection("feedstock_params.csv")
pt_options = csv2widget_collection("pretreatment_params.csv")
eh_options = csv2widget_collection("enzymatic_hydrolysis_params.csv")
br_options = csv2widget_collection("bioreactor_params.csv")


FS_model = Feedstock(fs_options)
PT_model = Pretreatment(pt_options, hpc_run)
EH_model = EnzymaticHydrolysis(eh_options, hpc_run)
BR_model = Bioreactor(br_options, hpc_run)


PT_model.run(verbose=False)

# EH_model.run_eh_cfd_surrogate()
EH_model.run_eh_cfd_simulation()
# EH_model.run_eh_lignocellulose_model()

# BR_model.run_biorector_cfd_surrogate()
BR_model.run_biorector_cfd_simulation()
