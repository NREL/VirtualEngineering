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
    display_name: Python 3 (ipykernel)
    language: python
    name: python3
---

---

# Virtual Engineering

The first step is to select "Cell" > "Run All" from the toolbar.  This will initialize all the widgets and allow you to interact with the unit operation options via the GUI controls.

<img src="docs/figures/three_unit_flow.png" alt="flowchart" width="800"/>

```python
from ipywidgets import *
from IPython.display import HTML, clear_output
import os
import sys
import numpy as np

#================================================================
# attempt to capture the parent directory in case of errors
if not 'notebookDir' in globals():
    notebookDir = os.getcwd()
#os.chdir(notebookDir)  # If you changed the current working dir, this will take you back to the workbook dir.
#================================================================

# imports from vebio modules
from vebio.WidgetFunctions import WidgetCollection, OptimizationWidget
from vebio.FileModifiers import write_file_with_replacements
from vebio.Utilities import get_host_computer, yaml_to_dict, dict_to_yaml
from vebio.RunFunctions import Pretreatment, Feedstock, EnzymaticHydrolysis, Bioreactor
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
fs_options = WidgetCollection()

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

initial_porosity_options = {'value': 0.8,
    'max': 1,
    'min': 0,
    'description': r'Initial Porosity',
    'description_tooltip': 'The initial porous fraction of the biomass particles.  Must be in the range [0, 1]'
}

fs_options.initial_porosity = OptimizationWidget('BoundedFloatText', initial_porosity_options)


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
pt_options = WidgetCollection()

initial_acid_conc_options = {'value': 0.0001,
    'max': 0.001,
    'min': 0.00005,
    'description': 'Acid Loading',
    'description_tooltip': 'The initial concentration of acid (mol/mL).  Must be in the range [0, 1]'}

#### this needs to be changed to g acid / g bone-dry biomass (then converted in the run function) ####
pt_options.initial_acid_conc = OptimizationWidget('BoundedFloatText', initial_acid_conc_options)

steam_temperature_options = {'value': 150.0,
    'max': 250.3,
    'min': 3.8,
    'description': 'Steam Temperature',
    'description_tooltip': r'The fixed temperature of the steam ($^\circ$C).',
}

pt_options.steam_temperature = OptimizationWidget('BoundedFloatText', steam_temperature_options)

# Conversion from celsius to kelvin
pt_options.steam_temperature.scaling_fn = lambda C : C + 273.15


initial_solid_fraction_options = {'value': 0.745,
    'max': 1,
    'min': 0,
    'description': r'Initial FIS$_0$',
    'description_tooltip': 'The initial fraction of insoluble solids (kg/kg).  Must be in the range [0, 1]'
}

pt_options.initial_solid_fraction = OptimizationWidget('BoundedFloatText', initial_solid_fraction_options)

final_time_options = {'value': 8.3,
    'max': 1440,
    'min': 1,
    'description': 'Final Time',
    'description_tooltip': r'Total simulation time (min).  Must be $\geq$ 1'
}

pt_options.final_time = OptimizationWidget('BoundedFloatText', final_time_options)

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
eh_options = WidgetCollection()

eh_options.model_type = widgets.RadioButtons(
    options = ['Lignocellulose Model', 'CFD Surrogate', 'CFD Simulation'],
    value = 'CFD Surrogate',
    description = 'Model Type',
    disabled = False,
    description_tooltip = 'Specifies the solver to use for the EH step, "CFD Simulation" requires HPC resources.'
)

lambda_e_options = {'value': 30.0,
    'max': 300.0,
    'min': 5.0,
    'description': 'Enzymatic Load',
    'description_tooltip': 'Ratio of the enzyme mass to the total solution mass (mg/g).  Must be in the range [0, 1000]'
}

eh_options.lambda_e = OptimizationWidget('BoundedFloatText', lambda_e_options)

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
br_options = WidgetCollection()

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

#================================================================

# Display the widgets
br_options.display_all_widgets()
```

---

### Choosing objective for optimization

```python
obj_widget = widgets.Dropdown(
    options=[('Biorector:            OUR',   ('bioreactor_outpit', 'our')), 
             ('Enzymatic Hydrolysis: rho_g', ('enzymatic_output', 'rho_g')), 
             ('Enzymatic Hydrolysis: rho_x', ('enzymatic_output', 'rho_x')),
             ('Enzymatic Hydrolysis: rho_sL',('enzymatic_output', 'rho_sL')),
             ('Enzymatic Hydrolysis: rho_f', ('enzymatic_output', 'rho_f')),
             ('Pretreatment:         fis_0', ('enzymatic_output', 'fis_0')),
             ('Pretreatment:         X_X',   ('enzymatic_output', 'X_X')),
             ('Pretreatment:         X_G',   ('enzymatic_output', 'X_G')),
             ('Pretreatment:         rho_x', ('enzymatic_output', 'rho_x')),
             ('Pretreatment:         rho_f', ('enzymatic_output', 'rho_f'))
            ],
    value=('bioreactor_outpit', 'our'),
    description='Optimization objective',
)



obj_widget = widgets.RadioButtons(
    options = ['our', 'rho_g', 'rho_x', 'rho_sL', 'rho_f', 'fis_0', 'X_X', 'X_G'],
    value = 'our',
    description = 'Optimization objective',
    disabled = False,
    description_tooltip = 'Specifies the objective to use in optimization.'
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
    
    with open('param_sweep.csv', 'w') as fp:
        fp.write('# Iteration, -, Enzyme Loading, OUR\n')
    
    sweep_ctr = 0
    nn = 4 # The number of points to select across each value
    
    initial_acid_conc_range = np.linspace(0.00005, 0.001, nn)   #[0.00005, 0.0001, 0.00015, 0.0002, 0.00025, 0.0003] 
#     initial_solid_fraction_range = [0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95]
#     final_time_range = [6, 7, 8, 9, 10, 11]
    enz_loading_range = np.linspace(5, 300, nn)    #[10, 20, 30, 40, 50, 60]
    
    # Set global paths and files for communication between operations
    os.chdir(notebookDir)
    params_filename = 'virteng_params_opt.yaml'
    
    FD_model = Feedstock(params_filename, fs_options)
    PT_model = Pretreatment(notebookDir, params_filename, pt_options)
    EH_model = EnzymaticHydrolysis(notebookDir, params_filename, eh_options, hpc_run)
    BR_model = Bioreactor(notebookDir, params_filename, br_options, hpc_run)
    
    for initial_acid_conc in initial_acid_conc_range:
        for enz_loading in enz_loading_range:
          # Set the swept value
            PT_model.initial_acid_conc=initial_acid_conc
            EH_model.lambda_e=enz_loading
            
            PT_model.run(verbose=False)         # Run the pretreatment model
            EH_model.run(verbose=False) # Run the enzymatic hydrolysis model
            BR_model.run(verbose=True)            # Run the bioreactor model

            sweep_ctr += 1

            with open('param_sweep.csv', 'a') as fp:
                v1 = PT_model.initial_acid_conc
#                 v2 = pt_options.initial_solid_fraction.widget.value
#                     v3 = pt_options.final_time.widget.value
                v3 = EH_model.lambda_e

                output_dict = yaml_to_dict(params_filename)
                obj = output_dict['bioreactor_output']['our']

#                 fp.write(f'{sweep_ctr:.0f}, {v1:.9e}, {v2:.9e}, {v3:.9e}, {obj:.9e}\n')
                fp.write(f'{sweep_ctr:.0f}, {v1:.9e}, {v3:.9e}, {obj:.9e}\n')

    
sweep_button.on_click(sweep_button_action)

#================================================================

```

```python
import matplotlib.pyplot as plt

param_sweep_fn = os.path.join(notebookDir, 'param_sweep.csv')
if os.path.exists(param_sweep_fn):
    sweeps = np.loadtxt(param_sweep_fn, delimiter=',', skiprows=1)
    extent = np.min(sweeps[:, 1]), np.max(sweeps[:, 1]), np.min(sweeps[:, 2]), np.max(sweeps[:, 2])
    nn = int(np.sqrt(len(sweeps)))
    OUR = sweeps[:, 3].reshape(nn, nn)

    shw = plt.imshow(OUR.T, extent=extent, aspect='auto', origin='lower')
    bar = plt.colorbar(shw)
    bar.set_label('OUR')
    plt.xlabel('initial_acid_conc')
    plt.ylabel('enz_loading')
```

 ## Optimize

Press the Optimize button below to launch the optimization of the start-to-finish operation using the above values as initial conditions.

This example **maximizes OUR** by **changing user-specified pretreatment options**.

```python
# import scipy.optimize as opt
from vebio.OptimizationFunctions import Optimization

#================================================================

opt_button = widgets.Button(
    description = 'Optimize.',
    tooltip = 'Optimize for OUR using the conditions above as an initial guess.',
    layout =  {'width': '200px', 'margin': '25px 0px 25px 170px'}, 
    button_style = 'warning'
)

#================================================================

# run_button_output = widgets.Output()
display(opt_button)

#================================================================
# Define a function to be executed each time the run button is pressed
def opt_button_action(b):
    clear_output()
    display(opt_button)
    
    params_filename = 'virteng_params_optimization.yaml'
    opt_results_file = 'optimization_results.csv'
    
    Opt = Optimization(fs_options, pt_options, eh_options, br_options, obj_widget,
                       hpc_run, notebookDir, params_filename, opt_results_file)
    
    opt_result = Opt.scipy_minimize(Opt.objective_function)
    print(opt_result)
    
opt_button.on_click(opt_button_action)

#================================================================
```

```python
opt_results = np.loadtxt(os.path.join(notebookDir, 'optimization_results.csv'), delimiter=',', skiprows=1)
print()
shw = plt.imshow(OUR.T, extent=extent, aspect='auto', origin='lower')
bar = plt.colorbar(shw)
bar.set_label('OUR')
plt.xlabel('Acid Loading')
plt.ylabel('Enzymatic Load')
plt.scatter(opt_results[:, 1], opt_results[:, 2], s=50, c='k', marker='o')
plt.plot(opt_results[:, 1], opt_results[:, 2], color='k')
plt.scatter(opt_results[-1, 1], opt_results[-1, 2], s=50, c='r', marker='o')
```

---

```python
a = HTML('''<script>
code_show=true; 
function code_toggle() {
 if (code_show){
 $('div.input').hide();
 } else {
 $('div.input').show();
 }
 code_show = !code_show
} 
$( document ).ready(code_toggle);
</script>
<form action="javascript:code_toggle()"><input type="submit" \
value="Toggle code visibility (hidden by default)."></form>''')

display(a)
```

```python
# # reload "run functions" code if needed
# if False:
#     from importlib import reload
#     import vebio.RunFunctions
#     reload(vebio.RunFunctions)
#     from vebio.RunFunctions import run_pretreatment, run_enzymatic_hydrolysis, run_bioreactor
```

```python

```
