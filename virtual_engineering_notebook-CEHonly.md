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
from vebio.RunFunctions import run_pretreatment, run_enzymatic_hydrolysis, run_bioreactor, run_CEH
# add path for no-CFD EH model
sys.path.append(os.path.join(notebookDir, "submodules/CEH_EmpiricalModel/src/core/"))


#================================================================
# See if we're running on Eagle or on a laptop
hpc_run = get_host_computer()
#================================================================
```

---

### 1. Continuous Enzymatic Hydrolysis Operation

Set the options for the continuous enzymatic hydrolysis operation using CEH model
#### Process flow chart
                                                makeup buffer1   makeup buffer 2
                                                      |                 |
                                                      v                 v
    DMR material -> Plug flow reactor EH stage -> CEH reactor 1 -> CEH reactor 2 -> exit stream

                                                      |                 |
                                                      v                 v
                                               sugar stream 1    sugar stream 2


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

ceh_options.fis_0 = widgets.BoundedFloatText(
    value = 0.20,
    max = 1.0,
    min = 0.0,
    description = r'Initial FIS$_0$',
    description_tooltip = 'The  initial fraction of insoluble solids of DMR material (kg/kg).  Must be in the range [0, 1]'
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

ceh_options.lignin_solid_fraction = widgets.BoundedFloatText(
    value = 0.113,
    max = 1.0,
    min = 0.0,
    description = 'Initial $X_L$',
    description_tooltip = 'The initial fraction of solids that is Xylan (kg/kg) after DMR pretreatment.  Must be in the range [0, 1]'
)

ceh_options.facile_fraction_glucan = widgets.BoundedFloatText(
    value = 0.6,
    max = 1.0,
    min = 0.0,
    description = 'Facile Glucan',
    description_tooltip = 'Facile glucan fraction.  Must be in the range [0, 1]'
)

ceh_options.t_final = widgets.BoundedFloatText(
    value = 25.0,
    min = 1.0,
    max = 100.0,
    description = 'Final Time',
    description_tooltip = r'The total time of the simulation (h).  Must be $\geq$ 1'
)

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
    #run_pretreatment(notebookDir, params_filename, fs_options, pt_options)
    
    # Run the enzymatic hydrolysis model
    #run_enzymatic_hydrolysis(notebookDir, params_filename, eh_options, hpc_run)
    
    # Run the CEH modelS
    run_CEH(notebookDir, params_filename, ceh_options)
    
    # Run the bioreactor model
    #run_bioreactor(notebookDir, params_filename, br_options, hpc_run)
    
run_button.on_click(run_button_action)

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
