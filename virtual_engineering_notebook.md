---
jupyter:
  jupytext:
    formats: ipynb,md
    text_representation:
      extension: .md
      format_name: markdown
      format_version: '1.2'
      jupytext_version: 1.4.0
  kernelspec:
    display_name: Python 3
    language: python
    name: python3
---

---

# Virtual Engineering

The first step is to select "Cell" > "Run All" from the toolbar.  This will initialize all the widgets and allow you to interact with the unit operation options via the GUI controls.

![flowchart](three_unit_flow.png)


```python
%reset -f
from ipywidgets import *
from IPython.display import HTML, clear_output
import yaml
import os
import time
from shutil import copyfile

import vebio.WidgetFunctions as wf
from vebio.FileModifiers import write_file_with_replacements
from vebio.Utilities import get_host_computer, yaml_to_dict, dict_to_yaml

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

# Create the collection of widgets
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
    description = r'initial porosity',
    description_tooltip = 'The initial forous fraction of the biomass particles.  Must be in the range [0, 1]'
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

# Create the collection of widgets
pt_options = wf.WidgetCollection()

pt_options.initial_acid_conc = widgets.BoundedFloatText(
    value = 0.0001,
    max = 1.0,
    min = 0.0,
    description = 'Acid Loading',
    description_tooltip = 'The initial concentration of acid (g/g).  Must be in the range [0, 1]'
)

pt_options.steam_temperature = widgets.BoundedFloatText(
    value = 423,
    max = 1000,
    min = 100,
    description = 'Steam Temperature',
    description_tooltip = r'The fixed temperature of the steam (K).'
)

pt_options.bulk_steam_conc = widgets.BoundedFloatText(
    value = 0.0001,
    max = 1.0,
    min = 0.0,
    description = 'Bulk Steam Concentration',
    description_tooltip = 'The ambient steam concentration.  Must be in the range [0, 1]'
)

pt_options.initial_solid_fraction = widgets.BoundedFloatText(
    value = 0.745,
    max = 1,
    min = 0,
    description = r'Initial FIS$_0$',
    description_tooltip = 'The initial fraction of insoluble solids (kg/kg).  Must be in the range [0, 1]'
)

pt_options.final_time = widgets.BoundedFloatText(
    value = 100,
    max = 1e16,
    min = 1,
    description = 'Final Time',
    description_tooltip = r'Total simulation time (s).  Must be $\geq$ 1'
)

pt_options.show_plots = widgets.Checkbox(
    value = False,
    description = 'Show Plots'
)


#================================================================

# Display the widgets
pt_options.display_all_widgets()

#================================================================

```

---

### 2. Enzymatic Hydrolysis Operation

Set the options for the enzymatic hydrolysis operation using a two-phase reaction rate model below.


```python
#================================================================

# Create the collection of widgets
eh_options = wf.WidgetCollection()

eh_options.lambda_e = widgets.BoundedFloatText(
    value = 0.03,
    max = 1.0,
    min = 0.0,
    description = 'Enzymatic Load',
    description_tooltip = 'Ratio of the enzyme mass to the total solution mass (kg/kg).  Must be in the range [0, 1]'
)

eh_options.fis_0_target = widgets.BoundedFloatText(
    value = 0.05,
    max = 1.0,
    min = 0.0,
    description = r'FIS$_0$ Target',
    description_tooltip = 'The target value for initial fraction of insoluble solids *after* dilution (kg/kg).  Must be in the range [0, 1]'
)

eh_options.t_final = widgets.BoundedFloatText(
    value = 100.0,
    min = 1.0,
    max = 1e16,
    description = 'Final Time',
    description_tooltip = r'The total time of the simulation (h).  Must be $\geq$ 1'
)

eh_options.show_plots = widgets.Checkbox(
    value = False,
    description = 'Show Plots'
)

#================================================================

# Display the widgets
eh_options.display_all_widgets()

#================================================================

```

---

### 3. Bioreaction Operation

Set the options for the bubble column bioreaction operation below.


```python
#================================================================

# Create the collection of widgets
br_options = wf.WidgetCollection()

br_options.t_final = widgets.BoundedFloatText(
    value = 100.0,
    min = 1.0,
    max = 1e16,
    description = 'Final Time',
    description_tooltip = r'The total time of the simulation (h).  Must be $\geq 1$'
                                    # is this really 'h'? current quasi-steady simulations only run 10s of seconds
)

#================================================================

# Display the widgets
br_options.display_all_widgets()

#================================================================

```

---

## Run Model

When finished setting options for all unit operations, press the button below to run the complete model.


```python
#================================================================

run_button = widgets.Button(
    description = 'Run All.',
    tooltip = 'Execute the model start-to-finish with the properties specified above.'
)

#================================================================

# run_button_output = widgets.Output()
display(run_button)

#================================================================

# Define a function to be executed each time the run button is pressed
def run_button_action(b):
    clear_output()
    display(run_button)
    
    # Store the current working directory
    parent_path = os.getcwd()
    
    # Run the pretreatment model
    # These print statements do not show for me until the run is complete. Something that should be fixed
    # at some point, JJS 3/22/20
    print('Running Pretreatment Model')
    os.chdir('pretreatment_model/test/')
    _ = fs_options.export_widgets_to_yaml('fs_input.yaml')
    _ = pt_options.export_widgets_to_yaml('pt_input.yaml', 'fs_input.yaml')

    %run ptrun.py 'pt_input.yaml' 'pt_output.yaml'

    if pt_options.show_plots.value:
        %run postprocess.py 'out_\*.dat' exptdata_150C_1acid.dat
    copyfile('pt_output.yaml', '../../two_phase_batch_model/%s' % ('pt_to_eh_input.yaml'))
    os.chdir(parent_path)
    print('\nFinished Pretreatment')
    
    # Run the enzymatic hydrolysis model
    print('\nRunning Enzymatic Hydrolysis Model')
    os.chdir('two_phase_batch_model/')
    _ = eh_options.export_widgets_to_yaml('eh_input.yaml', 'pt_to_eh_input.yaml')
    
    # output argument is not currently used, but it should be for passing info to the bioreaction 
    # unit operation, JJS 3/22/20
    %run two_phase_batch_model.py 'eh_input.yaml' 'eh_output.yaml'
    copyfile('eh_output.yaml', '../bioreactor/bubble_column/constant/%s' % ('eh_to_br_input.yaml'))
    os.chdir(parent_path)
    print('\nFinished Enzymatic Hydrolysis')
    
    # Run the bioreactor model
    os.chdir('bioreactor/bubble_column/constant/')

    #================================================================

    # FIXME: these lines of code should really exist in the fvOptions file itself
    # Create a dictionary of all the replacement values
    # where the key is a unique string to identify the definition
    # and the value is the corresponding value to assign

    # Make changes to the constant/fvOptions file
    
    eh_to_br_dict = yaml_to_dict('eh_to_br_input.yaml')        

    fvOptions_replacements = {}
    for k, v in eh_to_br_dict.items():
        fvOptions_replacements['double %s' % (k)] = v

    write_file_with_replacements('fvOptions_base', fvOptions_replacements)

    # Make changes to the system/controlDict file
    os.chdir('../system/')
    controlDict_replacements = {}
    controlDict_replacements['endTime '] = br_options.t_final.value
    
    write_file_with_replacements('controlDict_base', controlDict_replacements)
    
    #================================================================

    os.chdir('../')
    print('\nRunning Bioreactor Model')
    if hpc_run:
        # call function to update ovOptions # fvOptions?
        !sbatch ofoamjob
    else:
        print('Cannot run bioreactor without HPC resources.')
        print('$ sbatch ofoamjob')
        print(os.getcwd())
    print('\nFinished Bioreactor')

    # Return to the parent directory
    os.chdir(parent_path)
    
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
value="Toggle notebook code visibility (hidden by default)."></form>''')

display(a)
```

```python

```

```python

```
