---
jupyter:
  jupytext:
    formats: ipynb,md
    text_representation:
      extension: .md
      format_name: markdown
      format_version: '1.2'
      jupytext_version: 1.7.1
  kernelspec:
    display_name: Python 3
    language: python
    name: python3
---

---

# Virtual Engineering

The first step is to select "Cell" > "Run All" from the toolbar.  This will initialize all the widgets and allow you to interact with the unit operation options via the GUI controls.

<img src="three_unit_flow.png" alt="flowchart" width="800"/>

```python
from ipywidgets import *
from IPython.display import HTML, clear_output
import os
import sys

#================================================================
# attempt to capture the parent directory in case of errors
if not 'notebookDir' in globals():
    notebookDir = os.getcwd()
#os.chdir(notebookDir)  # If you changed the current working dir, this will take you back to the workbook dir.
#================================================================

# imports from vebio modules
import vebio.WidgetFunctions as wf
from vebio.FileModifiers import write_file_with_replacements
from vebio.Utilities import get_host_computer, yaml_to_dict, dict_to_yaml
from vebio.RunFunctions import run_pretreatment, run_enzymatic_hydrolysis, run_bioreactor
# add path for no-CFD EH model
sys.path.append(os.path.join(notebookDir, "submodules/CEH_EmpiricalModel/"))

#================================================================
# See if we're running on Eagle or on a laptop
hpc_run = get_host_computer()
#================================================================
```

<!-- #region -->
---

## Set Virtual Engineering Options


### 0. Feedstock properties

Set the feedstock properties.

<!-- #endregion -->

```python
#================================================================

# Create the collection of widgets for feedstock options
fs_options = wf.WidgetCollection()

fs_options.xylan_solid_fraction = widgets.BoundedFloatText(
    value = 0.263,
    max = 1,
    min = 0,
    description = r'Initial $X_X$',
    description_tooltip = 'The initial fraction of solids that is xylan (kg/kg).  Must be in the range [0, 1]'
)

fs_options.glucan_solid_fraction = widgets.BoundedFloatText(
    value = 0.40,
    max = 1,
    min = 0,
    description = r'Initial $X_G$',
    description_tooltip = 'The initial fraction of solids that is glucan (kg/kg).  Must be in the range [0, 1]'
)

fs_options.initial_porosity = widgets.BoundedFloatText(
    value = 0.8,
    max = 1,
    min = 0,
    description = r'Initial Porosity',
    description_tooltip = 'The initial porous fraction of the biomass particles.  Must be in the range [0, 1]'
)


#================================================================

# Display the widgets
fs_options.display_all_widgets()

#================================================================

```

---

### 1. Pretreatment Operation

Set the options for the pretreatment operation below.

```python
#================================================================

# Create the collection of widgets for pretreatment options
pt_options = wf.WidgetCollection()

#### this needs to be changed to g acid / g bone-dry biomass (then converted in the run function) ####
pt_options.initial_acid_conc = widgets.BoundedFloatText(
    value = 0.0001,
    max = 1.0,
    min = 0.0,
    description = 'Acid Loading',
    description_tooltip = 'The initial concentration of acid (mol/mL).  Must be in the range [0, 1]'
)

pt_options.steam_temperature = widgets.BoundedFloatText(
    value = 150.0,
    max = 250.3,
    min = 3.8,
    description = 'Steam Temperature',
    description_tooltip = r'The fixed temperature of the steam ($^\circ$ C).',
)
# Conversion from celsius to kelvin
pt_options.steam_temperature.scaling_fn = lambda C : C + 273.15

# This has been replaced with a value from a saturated steam
# lookup table based on temperature in the pretreatment run function.
# pt_options.bulk_steam_conc = widgets.BoundedFloatText(
#     value = 0.0001,
#     max = 1.0,
#     min = 0.0,
#     description = 'Bulk Steam Concentration',
#     description_tooltip = 'The ambient steam concentration (mol/mL).  Must be in the range [0, 1]'
# )

pt_options.initial_solid_fraction = widgets.BoundedFloatText(
    value = 0.745,
    max = 1,
    min = 0,
    description = r'Initial FIS$_0$',
    description_tooltip = 'The initial fraction of insoluble solids (kg/kg).  Must be in the range [0, 1]'
)

pt_options.final_time = widgets.BoundedFloatText(
    value = 8.3, #100, #2400
    max = 1440,
    min = 1,
    description = 'Final Time',
    description_tooltip = r'Total simulation time (min).  Must be $\geq$ 1'
)
# Conversion from minutes to seconds
pt_options.final_time.scaling_fn = lambda s : 60.0 * s

pt_options.show_plots = widgets.Checkbox(
    value = False,
    description_tooltip = 'Show Plots'
)


#================================================================

# Display the widgets
pt_options.display_all_widgets()

#================================================================

```

---

### 2. Enzymatic Hydrolysis Operation

Set the options for the enzymatic hydrolysis operation using either a two-phase reaction rate model or high-fidelity CFD below.


```python
#================================================================

# Create the collection of widgets
eh_options = wf.WidgetCollection()

eh_options.model_type = widgets.RadioButtons(
    options = ['Lignocellulose Model', 'CFD Surrogate', 'CFD Simulation'],
    value = 'CFD Surrogate',
    description = 'Model Type',
    disabled = False,
    description_tooltip = 'Specifies the solver to use for the EH step, "CFD Simulation" requires HPC resources.'
)

eh_options.lambda_e = widgets.BoundedFloatText(
    value = 30.0,
    max = 1000.0,
    min = 0.0,
    description = 'Enzymatic Load',
    description_tooltip = 'Ratio of the enzyme mass to the total solution mass (mg/g).  Must be in the range [0, 1000.0]'
)
# Conversion from mg/g to kg/kg
eh_options.lambda_e.scaling_fn = lambda e : 0.001 * e

eh_options.fis_0 = widgets.BoundedFloatText(
    value = 0.05,
    max = 1.0,
    min = 0.0,
    description = r'FIS$_0$ Target',
    description_tooltip = 'The target value for initial fraction of insoluble solids *after* dilution (kg/kg).  Must be in the range [0, 1]'
)

eh_options.t_final = widgets.BoundedFloatText(
    value = 24.0,
    min = 1.0,
    max = 24.0,
    description = 'Final Time',
    description_tooltip = r'The total time of the simulation (h).  Must be $\geq$ 1'
)

eh_options.show_plots = widgets.Checkbox(
    value = False,
    description_tooltip = 'Show Plots',
    disabled = True,
)

#================================================================

# Display the widgets
eh_options.display_all_widgets()

#================================================================

def model_type_action(change):
        
    if eh_options.model_type.value == 'Lignocellulose Model':
        # Lignocellulose Model
        eh_options.show_plots.value = False
        eh_options.show_plots.disabled = False
        eh_options.show_plots.description_tooltip = 'Show Plots'

    else:
        # Surrogate Model, CFD Simulation
        eh_options.show_plots.value = False
        eh_options.show_plots.disabled = True
        eh_options.show_plots.description_tooltip = 'Show Plots (Only available for lignocellulose model)'

eh_options.model_type.observe(model_type_action, names='value')
```

---

### 3. Bioreaction Operation

Set the options for the bubble column bioreaction operation below.


```python
#================================================================

# Create the collection of widgets
br_options = wf.WidgetCollection()

br_options.model_type = widgets.RadioButtons(
    options = ['CFD Surrogate', 'CFD Simulation'],
    value = 'CFD Surrogate',
    description = 'Model Type',
    disabled = False,
    description_tooltip = 'Specifies the solver to use for the bioreaction step, "CFD Simulation" requires HPC resources.'
)

br_options.t_final = widgets.BoundedFloatText(
    value = 100.0, # default 500
    min = 1.0,
    max = 1e16,
    description = 'Final Time',
    description_tooltip = r'The total time of the simulation (h).  Must be $\geq 1$'
                                    # is this really 'h'? current quasi-steady simulations only run tens of seconds
)

# br_options.use_cfd = widgets.Checkbox(
#     value = False,
#     description_tooltip = 'Use High-Fidelity CFD (Requires HPC Resources)',
# )

#================================================================

# Display the widgets
br_options.display_all_widgets()
```

---

## Run Model

When finished setting options for all unit operations, press the button below to run the complete model.


```python
#================================================================

run_button = widgets.Button(
    description = 'Run All.',
    tooltip = 'Execute the model start-to-finish with the properties specified above.',
    layout =  {'width': '200px', 'margin': '25px 0px 100px 170px'}, 
    button_style = 'success'
)

#================================================================

# run_button_output = widgets.Output()
display(run_button)

#================================================================

# Define a function to be executed each time the run button is pressed
def run_button_action(b):
    clear_output()
    display(run_button)
    
    # Set global paths and files for communication between operations
    os.chdir(notebookDir)
    params_filename = 'virteng_params.yaml'
    
    # Run the pretreatment model
    run_pretreatment(notebookDir, params_filename, fs_options, pt_options)
    
    # Run the enzymatic hydrolysis model
    run_enzymatic_hydrolysis(notebookDir, params_filename, eh_options, hpc_run)
    
    # Run the bioreactor model
    run_bioreactor(notebookDir, params_filename, br_options, hpc_run)
    
run_button.on_click(run_button_action)

#================================================================

```

---

```python
# INPUTS
# Gas velocity = 0.08
# Column height = 40.00
# Column diameter = 5.00
# Bubble diameter = 0.0060
# OUR_max = 88.13
# t_final = 100.0
# FINAL OUTPUTS (at t = 100.0 seconds)
# OUR = 0.0575

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
# reload "run functions" code if needed
if False:
    from importlib import reload
    import vebio.RunFunctions
    reload(vebio.RunFunctions)
    from vebio.RunFunctions import run_pretreatment, run_enzymatic_hydrolysis, run_bioreactor
```

```python
def temp(a, parent=None):
    print(a, parent)
    
    
temp(4, parent=6)
```

```python
import numpy as np
from scipy import interpolate as interp
steam = np.genfromtxt('sat_steam_table.csv', delimiter=',', skip_header=1)
# print(np.shape(steam))
# print(steam)

x = interp.interp1d(steam[:, 1], steam[:, 3])
y = x(150)
print(y, y/18.01528/1000.0)
```

```python

```
