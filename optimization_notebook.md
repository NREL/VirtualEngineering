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
    params_filename = 'virteng_params_opt.yaml'
    # Run the pretreatment model
    run_pretreatment(notebookDir, params_filename, fs_options, pt_options)
    
    # Run the enzymatic hydrolysis model
    run_enzymatic_hydrolysis(notebookDir, params_filename, eh_options, hpc_run)
    
    # Run the bioreactor model
    run_bioreactor(notebookDir, params_filename, br_options, hpc_run)
    
run_button.on_click(run_button_action)

#================================================================

```

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
#         fp.write('# Iteration, Acid Loading, Initial FIS_0, Final Time, OUR\n')
        fp.write('# Iteration, Acid Loading, Enzyme Loading, OUR\n')
    
    sweep_ctr = 0
    
    nn = 3 # The number of points to select across each value
    
    initial_acid_conc_range = np.linspace(0.00005, 0.001, 12)#[0.00005, 0.0001, 0.00015, 0.0002, 0.00025, 0.0003]
    
#     initial_solid_fraction_range = [0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95]

#     final_time_range = [6, 7, 8, 9, 10, 11]
    enz_loading_range = np.linspace(5, 300, 12)#[10, 20, 30, 40, 50, 60]
    
    for initial_acid_conc in initial_acid_conc_range:
        for enz_loading in enz_loading_range:
            # Set global paths and files for communication between operations
            os.chdir(notebookDir)
            params_filename = 'virteng_params_opt.yaml'

            # Set the swept value
            pt_options.initial_acid_conc.widget.value = initial_acid_conc
#             pt_options.initial_solid_fraction.widget.value = initial_solid_fraction
#                 pt_options.final_time.widget.value = final_time
            eh_options.lambda_e.widget.value = enz_loading

            # Run the pretreatment model
            run_pretreatment(notebookDir, params_filename, fs_options, pt_options, verbose=False)

            # Run the enzymatic hydrolysis model
            run_enzymatic_hydrolysis(notebookDir, params_filename, eh_options, hpc_run, verbose=False)

            # Run the bioreactor model
            run_bioreactor(notebookDir, params_filename, br_options, hpc_run, verbose=True)

            sweep_ctr += 1

            with open('param_sweep.csv', 'a') as fp:
                v1 = pt_options.initial_acid_conc.widget.value
#                 v2 = pt_options.initial_solid_fraction.widget.value
#                     v3 = pt_options.final_time.widget.value
                v3 = eh_options.lambda_e.widget.value

                output_dict = yaml_to_dict(params_filename)
                obj = output_dict['bioreactor_output']['our']

#                 fp.write(f'{sweep_ctr:.0f}, {v1:.9e}, {v2:.9e}, {v3:.9e}, {obj:.9e}\n')
                fp.write(f'{sweep_ctr:.0f}, {v1:.9e}, {v3:.9e}, {obj:.9e}\n')

    
sweep_button.on_click(sweep_button_action)

#================================================================

```

 ## Optimize

Press the Optimize button below to launch the optimization of the start-to-finish operation using the above values as initial conditions.

This example **maximizes OUR** by **changing user-specified pretreatment options**.

```python
import scipy.optimize as opt

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

def opt_callback(free_variables):
    pass
#     print('Controls:', free_variables)

fn_evals = 0
objective_scaling = 1.0

def objective_function(free_variables, notebookDir, params_filename,
                 fs_options, pt_options, eh_options, br_options, hpc_run, opt_results_file, x_0_names):
            
    global fn_evals
    global objective_scaling
    
    ctr = 0
    
    dimensional_values = []
    
    # Update the controls with the latest values
    for wc in [fs_options, pt_options, eh_options]:
        for widget_name, widget in wc.__dict__.items():        
            if isinstance(widget, OptimizationWidget) and widget.is_control.value == True:
                lb = widget.widget.min
                ub = widget.widget.max

                dimensional_value = free_variables[ctr]*(ub-lb) + lb
                dimensional_values.append(dimensional_value)

                widget.widget.value = dimensional_value
                ctr += 1
    
#     for item in fs_options.__dict__.items():
#         w_name = item[0]
#         w = item[1]
        
#         if hasattr(w, 'optimize'):
#             if getattr(w, 'optimize') == True:
#                 setattr(w, 'value', free_variables[ctr])
#                 ctr += 1
                
    # Set global paths and files for communication between operations
    os.chdir(notebookDir)
    
    # Turn off printed outputs from unit operations
    if fn_evals == 0:
        v_flag = True
    else:
        v_flag = False
    
    # Run the pretreatment model
    run_pretreatment(notebookDir, params_filename, fs_options, pt_options, verbose=v_flag)
    
    # Run the enzymatic hydrolysis model
    run_enzymatic_hydrolysis(notebookDir, params_filename, eh_options, hpc_run, verbose=v_flag)
    
    # Run the bioreactor model
    run_bioreactor(notebookDir, params_filename, br_options, hpc_run, verbose=v_flag)
    
    # Read the outputs into a dictionary
    output_dict = yaml_to_dict(params_filename)
    
    # The objective function in this case is OUR, we take the negative
    # so the minimize function sees the correct orientation
    obj = -output_dict['bioreactor_output']['our']        

    if fn_evals == 0:
        print('Beginning Optimization')
        objective_scaling = -1.0/obj
        
    fn_evals += 1
        
    with open(opt_results_file, 'a') as fp:
        fp.write('%d, ' % (fn_evals))
        for dv in dimensional_values:
            fp.write('%.15e, ' % (dv))
        fp.write('%.15e\n' % (obj))

    print('Iter = %3d: ' % (fn_evals), end='')
    for k, dv in enumerate(dimensional_values):
        print('%s = %12.9e, ' % (x_0_names[k], dv), end='')
    print('Objective = %12.9e' % (obj))
    
    obj *= objective_scaling
    print(obj)
    
    return obj


# Define a function to be executed each time the run button is pressed
def opt_button_action(b):
    clear_output()
    display(opt_button)
    
    assert br_options.model_type.value == 'CFD Surrogate'
    assert eh_options.model_type.value == 'CFD Surrogate'

    x_0 = []
    x_0_names = []
    x_0_bounds = []
    
    for wc in [fs_options, pt_options, eh_options]:
        for widget_name, widget in wc.__dict__.items():        
            if isinstance(widget, OptimizationWidget) and widget.is_control.value == True:
                print('Optimizing %s.' % widget_name)
                
                lb = widget.widget.min
                ub = widget.widget.max

                scaled_value = (widget.widget.value - lb)/(ub-lb)
                
#                 current_val = widget.widget.value

                # Here, we use the values of the controls scaled to the range [0, 1]
                x_0.append(scaled_value)
                x_0_names.append(widget_name)
#                 x_0_bounds.append((lb, ub))
                x_0_bounds.append((0.0, 1.0))
    
    if len(x_0) == 0:
        raise ValueError('No controls have been specified, retry with >= 1 control variables.')
    
    params_filename = 'virteng_params_optimization.yaml'
    
    
    opt_results_file = 'optimization_results.csv'
    with open(opt_results_file, 'w') as fp:
        fp.write('# Iteration, ')
        for control in x_0_names:
            fp.write('%s, ' % (control))
            
        fp.write('Objective\n')

     #L-BFGS-B
    global fn_evals
    global objective_scaling
    fn_evals = 0
    objective_scaling=1.0
    
    opt_result = opt.minimize(objective_function,
                 x_0,
                 (notebookDir, params_filename,
                  fs_options, pt_options,
                  eh_options, br_options,
                  hpc_run, opt_results_file, x_0_names),
                 method='SLSQP',
                 bounds=x_0_bounds,
                 callback=opt_callback)
        
    print(opt_result)
    
opt_button.on_click(opt_button_action)

#================================================================
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
