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
import vebio.WidgetFunctions as wf
from vebio.FileModifiers import write_file_with_replacements
from vebio.Utilities import get_host_computer, yaml_to_dict, dict_to_yaml
from vebio.RunFunctions import run_pretreatment, run_enzymatic_hydrolysis, run_bioreactor, run_CEH, run_ve_tea
# add path for no-CFD EH model
sys.path.append(os.path.join(notebookDir, "submodules/CEH_EmpiricalModel/src/core/"))
sys.path.append(os.path.join(notebookDir, "submodules/Aspen_tool/AutoAspen/"))


#================================================================
# See if we're running on Eagle or on a laptop
hpc_run = get_host_computer()
#================================================================
```

---

### 1. Continuous Enzymatic Hydrolysis Operation

Set the options for the continuous enzymatic hydrolysis operation using CEH model
#### Process flow chart
                    makeup buffer 1   makeup buffer 2  makeup buffer 3
                          |                 |                |
                          v                 v                v
    DMR material ->  CEH reactor 1 -> CEH reactor 2 -> CEH reactor 3 -> exit stream

                          |                 |                |
                          v                 v                v
                    sugar stream 1    sugar stream 2   sugar stream 3


```python
#================================================================

# Create the collection of widgets
ceh_options = wf.WidgetCollection()

ceh_options.model_type = widgets.RadioButtons(
    options = ['CEH Model'],
    value = 'CEH Model',
    description = 'Model Type',
    disabled = False,
    description_tooltip = 'Only have one option for now'
)

ceh_options.lambda_e = widgets.BoundedFloatText(
    value = 12.0,
    max = 1000.0,
    min = 0.0,
    description = 'Enzyme Loading',
    description_tooltip = 'Ratio of the enzyme mass to the total solution mass (mg/g).  Must be in the range [0, 1000]'
)
# Conversion from mg/g to kg/kg
ceh_options.lambda_e.scaling_fn = lambda e : 0.001 * e

ceh_options.f1_is = widgets.BoundedFloatText(
    value = 0.10,
    max = 1.0,
    min = 0.0,
    description = r'Inflow f1_is',
    description_tooltip = 'Inflow: The  initial fraction of insoluble solids in input stream for DMR slurry (kg/kg).  Must be in the range [0, 1]'
)

ceh_options.inflow_mass_flowrate = widgets.BoundedFloatText(
    value = 100.0,
    max = 1.0e8,
    min = 0.0,
    description = r'Inflow F1',
    description_tooltip = 'Inflow: mass flow rate for DMR slurry (kg/h).  Must be in the range [0, 1e8]'
)

ceh_options.glucan_solid_fraction = widgets.BoundedFloatText(
    value = 0.513,
    max = 1.0,
    min = 0.0,
    description = r'Initial $X_G$',
    description_tooltip = 'The initial fraction of solids that is Glucan (kg/kg) after DMR pretreatment.  Must be in the range [0, 1]'
)

ceh_options.xylan_solid_fraction = widgets.BoundedFloatText(
    value = 0.2535,
    max = 1.0,
    min = 0.0,
    description = r'Initial $X_X$',
    description_tooltip = 'The initial fraction of solids that is Xylan (kg/kg) after DMR pretreatment.  Must be in the range [0, 1]'
)

ceh_options.residue_soluble_lignin = widgets.BoundedFloatText(
    value = 2.0,
    max = 1000.0,
    min = 0.0,
    description = r'Residue $C_L$',
    description_tooltip = 'The residule soluble lignin in material inflow (g/L) of DMR slurry.  Must be in the range [0, 1000]'
)

#ceh_options.lignin_solid_fraction = widgets.BoundedFloatText(
#    value = 0.113,
#    max = 1.0,
#    min = 0.0,
#    description = 'Initial $X_L$',
#    description_tooltip = 'The initial fraction of solids that is Xylan (kg/kg) after DMR pretreatment.  Must be in the range [0, 1]'
#)

ceh_options.facile_fraction_glucan = widgets.BoundedFloatText(
    value = 0.6,
    max = 1.0,
    min = 0.0,
    description = 'Facile Glucan',
    description_tooltip = 'Facile glucan fraction.  Must be in the range [0, 1]'
)

ceh_options.fis_1 = widgets.BoundedFloatText(
    value = 0.07,
    max = 1.0,
    min = 0.0,
    description = r'Target FIS$_1$',
    description_tooltip = 'The target insoluble solids in CEH reactor 1 (kg/kg).  Must be in the range [0, 1]'
)

ceh_options.fis_2 = widgets.BoundedFloatText(
    value = 0.07,
    max = 1.0,
    min = 0.0,
    description = r'Target FIS$_2$',
    description_tooltip = 'The target insoluble solids in CEH reactor 2 (kg/kg).  Must be in the range [0, 1]'
)

ceh_options.fis_3 = widgets.BoundedFloatText(
    value = 0.07,
    max = 1.0,
    min = 0.0,
    description = r'Target FIS$_3$',
    description_tooltip = 'The target insoluble solids in CEH reactor 3 (kg/kg).  Must be in the range [0, 1]'
)

ceh_options.target_conv_1 = widgets.BoundedFloatText(
    value = 0.50,
    max = 1.0,
    min = 0.0,
    description = r'Target x$_1$',
    description_tooltip = 'The target hydrocarbon conversion in CEH reactor 1.  Must be in the range [0, 1]'
)

ceh_options.target_conv_2 = widgets.BoundedFloatText(
    value = 0.60,
    max = 1.0,
    min = 0.0,
    description = r'Target x$_2$',
    description_tooltip = 'The target hydrocarbon conversion in CEH reactor 2.  Must be in the range [0, 1]'
)

ceh_options.target_conv_3 = widgets.BoundedFloatText(
    value = 0.70,
    max = 1.0,
    min = 0.0,
    description = r'Target x$_3$',
    description_tooltip = 'The target hydrocarbon conversion in CEH reactor 3.  Must be in the range [0, 1]'
)

ceh_options.theta2_1 = widgets.BoundedFloatText(
    value = 0.50,
    max = 1.0,
    min = 0.0,
    description = r'Target theta2$_1$',
    description_tooltip = 'The target theta2 (makeup/feed ratio) in CEH reactor 1.  Must be in the range [0, 1]'
)

ceh_options.theta2_2 = widgets.BoundedFloatText(
    value = 0.50,
    max = 1.0,
    min = 0.0,
    description = r'Target theta2$_2$',
    description_tooltip = 'The target theta2 (makeup/feed ratio) in CEH reactor 2.  Must be in the range [0, 1]'
)

ceh_options.theta2_3 = widgets.BoundedFloatText(
    value = 0.50,
    max = 1.0,
    min = 0.0,
    description = r'Target theta2$_3$',
    description_tooltip = 'The target theta2 (makeup/feed ratio) in CEH reactor 3.  Must be in the range [0, 1]'
)



#ceh_options.t_final = widgets.BoundedFloatText(
#    value = 25.0,
#    min = 1.0,
#    max = 100.0,
#    description = 'Final Time',
#    description_tooltip = r'The total time of the simulation (h).  Must be $\geq$ 1'
#)

ceh_options.show_plots = widgets.Checkbox(
    value = False,
    description_tooltip = 'Show Plots',
)


#================================================================

# Display the widgets
ceh_options.display_all_widgets()

#================================================================

```

---

### 2. Techno-Economic Analysis Options


```python

```

```python
#================================================================

# Create the collection of widgets
tea_options = wf.WidgetCollection()

tea_options.aspen_filename = widgets.Text(
    value = 'bc1707a-sugars_CEH.bkp',
    description = 'Aspen File',
    disabled = True,
    description_tooltip = 'Path to Aspen backup file.'
)

tea_options.excel_filename = widgets.Text(
    value = 've_bc1707a-sugars_CEH.xlsm',
    description = 'Excel File',
    disabled = True,
    description_tooltip = 'Path to Excel calculation file.'
)


#================================================================

# Display the widgets
tea_options.display_all_widgets()

#================================================================

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
def run_button_action(b, clear_output=True):
    if clear_output:
        clear_output()
        display(run_button)

    # Set global paths and files for communication between operations
    os.chdir(notebookDir)
    params_filename = 'virteng_params.yaml'

    # Run the pretreatment model
    #run_pretreatment(notebookDir, params_filename, fs_options, pt_options)

    # Run the enzymatic hydrolysis model
    #run_enzymatic_hydrolysis(notebookDir, params_filename, eh_options, hpc_run)

    # Run the CEH modelS
    aspen_flowrate = 282628.0 # (kg/hr) for an FIS of 0.198
    aspen_fis = 0.198
    ceh_options.inflow_mass_flowrate.value = aspen_flowrate*aspen_fis/ceh_options.f1_is.value

    run_CEH(notebookDir, params_filename, ceh_options)

    # Run the bioreactor model
    #run_bioreactor(notebookDir, params_filename, br_options, hpc_run)
    
    run_ve_tea(notebookDir, params_filename, tea_options, verbose=True)

run_button.on_click(run_button_action)

#================================================================

#================================================================

run_sweep_button = widgets.Button(
    description = 'Run Param Sweep of FIS.',
    tooltip = 'Execute the model in parameter sweep mode with the properties specified above.',
    layout =  {'width': '200px', 'margin': '25px 0px 100px 170px'},
    button_style = 'success'
)

#================================================================

# run_button_output = widgets.Output()
display(run_sweep_button)

#================================================================

# Define a function to be executed each time the run button is pressed
def run_sweep_button_action(b):

    swept_param = ceh_options.f1_is

    nn = 10
    min_swept_param = 0.01
    max_swept_param = 0.1

    swept_vals = np.linspace(min_swept_param, max_swept_param, nn)

    for case_num, val in enumerate(swept_vals):
        print('================================================================')
        print(f'Sweep {case_num+1} of {nn}, {swept_param.description} = {val} ')
        print('================================================================')
        swept_param.value = val

        run_button_action(b, clear_output=False)


run_sweep_button.on_click(run_sweep_button_action)

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
