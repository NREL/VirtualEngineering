---
jupyter:
  jupytext:
    formats: ipynb,md
    text_representation:
      extension: .md
      format_name: markdown
      format_version: '1.3'
      jupytext_version: 1.14.4
  kernelspec:
    display_name: Python 3 (ipykernel)
    language: python
    name: python3
---

# Virtual Engineering with Optimization

The first step is to select "Cell" > "Run All" from the toolbar.  This will initialize all the widgets and allow you to interact with the unit operation options via the GUI controls.

```python
from IPython.display import Image
import os
Image(os.path.join(os.path.dirname("__file__"), 'docs', 'figures', 'three_unit_flow.png'), width=800)
```

```python
from ipywidgets import *
from IPython.display import HTML, clear_output
import os
import numpy as np

# imports from vebio modules
from vebio.WidgetFunctions import WidgetCollection, OptimizationWidget, csv2widget_collection
from vebio.Utilities import get_host_computer
from vebio.RunFunctions import Pretreatment, Feedstock, EnzymaticHydrolysis, Bioreactor
from vebio.OptimizationFunctions import Optimization
#================================================================
# See if we're running on Eagle or on a laptop
hpc_run = get_host_computer()
#================================================================
import logging
logger = logging.getLogger()
logger.setLevel(logging.CRITICAL)
```

## Set Virtual Engineering Options
### 0. Feedstock properties
Set the feedstock properties.

```python
fs_options = csv2widget_collection("feedstock_params.csv")
fs_options.display_all_widgets()
```

### 1. Pretreatment Operation

Set the options for the pretreatment operation below.

```python
pt_options = csv2widget_collection("pretreatment_params.csv")
pt_options.display_all_widgets()
```

---

### 2. Enzymatic Hydrolysis Operation

Set the options for the enzymatic hydrolysis operation using either a two-phase reaction rate model or high-fidelity CFD below.


```python
eh_options = csv2widget_collection("enzymatic_hydrolysis_params.csv")
eh_options.display_all_widgets()
```

---

### 3. Bioreaction Operation

Set the options for the bubble column bioreaction operation below.


```python
br_options = csv2widget_collection("bioreactor_params.csv")
br_options.display_all_widgets()
```

---

### Choosing objective for optimization

```python
obj_widget = widgets.Dropdown(
    options=[('Biorector:            OUR',   ('br_out', 'our')), 
             ('Enzymatic Hydrolysis: rho_g', ('eh_out', 'rho_g')), 
             ('Enzymatic Hydrolysis: rho_x', ('eh_out', 'rho_x')),
             ('Enzymatic Hydrolysis: rho_sL',('eh_out', 'rho_sL')),
             ('Enzymatic Hydrolysis: rho_f', ('eh_out', 'rho_f')),
             ('Pretreatment:         fis_0', ('eh_out', 'fis_0')),
             ('Pretreatment:         X_X',   ('eh_out', 'X_X')),
             ('Pretreatment:         X_G',   ('eh_out', 'X_G')),
             ('Pretreatment:         rho_x', ('eh_out', 'rho_x')),
             ('Pretreatment:         rho_f', ('eh_out', 'rho_f'))
            ],
    value=('br_out', 'our'),
    description='Objective:',
    tooltip = 'Specifies the objective to use in optimization.'
)

display(obj_widget)
```

---

### Parameter Sweeps

```python
#================================================================
sweep_button = widgets.Button(
    description = 'Run Sweep.',
    tooltip = 'Execute the model start-to-finish with the properties specified above.',
    layout =  {'width': '200px', 'margin': '25px 0px 100px 170px'}, 
    button_style = 'primary'
)
#================================================================
# run_button_output = widgets.Output()
display(sweep_button)
#================================================================

# Define a function to be executed each time the run button is pressed
def sweep_button_action(b):
    clear_output()
    display(sweep_button)
    Opt = Optimization(fs_options, pt_options, eh_options, br_options, obj_widget, hpc_run)
    Opt.parameter_grid_sweep(nn=10, results_file='sweep_params.csv')
    
sweep_button.on_click(sweep_button_action)
#===============================================================
```

```python
import matplotlib.pyplot as plt

param_sweep_fn = 'sweep_params.csv'
if os.path.exists(param_sweep_fn):
    with open(param_sweep_fn, 'r') as f:
        firstline = f.readline().split(',')
    Opt = Optimization(fs_options, pt_options, eh_options, br_options, obj_widget, hpc_run)
    sweeps = np.loadtxt(param_sweep_fn, delimiter=',', skiprows=1)
    bounds = Opt.var_real_bounds
    extent = bounds[0][0], bounds[0][1], bounds[1][0], bounds[1][1]
    nn = int(np.sqrt(len(sweeps)))
    OUR = sweeps[:, 3].reshape(nn, nn)
    shw = plt.imshow(OUR.T, extent=extent, aspect='auto', origin='lower')
    _cs2 = plt.contour(OUR.T, levels=np.arange(65, 69, 0.5), extent=extent, origin='lower', colors='blue')
    bar = plt.colorbar(shw)
    bar.add_lines(_cs2)
    bar.set_label(firstline[-1])
    plt.xlabel(firstline[1])
    plt.ylabel(firstline[2])
```

 ## Optimize

Press the Optimize button below to launch the optimization of the start-to-finish operation using the above values as initial conditions.

This example **maximizes OUR** by **changing user-specified pretreatment options**.

```python
logger.setLevel(logging.CRITICAL)

opt_button = widgets.Button(
    description = 'Optimize.',
    tooltip = 'Optimize for OUR using the conditions above as an initial guess.',
    layout =  {'width': '200px', 'margin': '25px 0px 25px 170px'}, 
    button_style = 'warning'
)
#================================================================

display(opt_button)

#================================================================
# Define a function to be executed each time the run button is pressed
def opt_button_action(b):
    clear_output()
    display(opt_button)
    
    params_filename = 'virteng_params_optimization.yaml'
    opt_results_file = 'optimization_results.csv'
    
    Opt = Optimization(fs_options, pt_options, eh_options, br_options, obj_widget,
                       hpc_run)
    
    opt_result = Opt.scipy_minimize(Opt.objective_function, opt_results_file=opt_results_file)
    print(opt_result)
    
opt_button.on_click(opt_button_action)
#================================================================
```

```python
opt_results = np.loadtxt('optimization_results.csv', delimiter=',', skiprows=1)
shw = plt.imshow(OUR.T, extent=extent, aspect='auto', origin='lower')
bar = plt.colorbar(shw)
bar.set_label('OUR')
plt.xlabel(firstline[1])
plt.ylabel(firstline[2])
plt.scatter(opt_results[:, 1], opt_results[:, 2], s=50, c='k', marker='o')
plt.plot(opt_results[:, 1], opt_results[:, 2], color='k')
plt.scatter(opt_results[-1, 1], opt_results[-1, 2], s=50, c='r', marker='o')
```

---

```python
# a = HTML('''<script>
# code_show=true; 
# function code_toggle() {
#  if (code_show){
#  $('div.input').hide();
#  } else {
#  $('div.input').show();
#  }
#  code_show = !code_show
# } 
# $( document ).ready(code_toggle);
# </script>
# <form action="javascript:code_toggle()"><input type="submit" \
# value="Toggle code visibility (hidden by default)."></form>''')

# display(a)
```

```python
# # reload "run functions" code if needed
# if False:
#     from importlib import reload
#     import vebio.RunFunctions
#     reload(vebio.RunFunctions)
#     from vebio.RunFunctions import run_pretreatment, run_enzymatic_hydrolysis, run_bioreactor
```
