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

eh_options.fis_0 = widgets.BoundedFloatText(
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

eh_options.use_cfd = widgets.Checkbox(
    value = False,
    description = 'Use High-Fidelity CFD (Requires HPC Resources)',
)

#================================================================

# Display the widgets
eh_options.display_all_widgets()

#================================================================

def use_cfd_action(b):
    if eh_options.use_cfd.value:
        eh_options.show_plots.value = False
        eh_options.show_plots.disabled = True
    else:
        eh_options.show_plots.value = False
        eh_options.show_plots.disabled = False
        

eh_options.use_cfd.observe(use_cfd_action)
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
    tooltip = 'Execute the model start-to-finish with the properties specified above.',
    layout =  {'width': '200px', 'margin': '25px 0px 100px 170px'}, 
    button_style = 'success'
)

#================================================================

# run_button_output = widgets.Output()
display(run_button)

#================================================================

def run_pretreatment():
    print('Running Pretreatment Model')
    parent_path = os.getcwd()
    
    # Export the feedstock and pretreatment options to a global yaml file
    fs_options.export_widgets_to_yaml('feedstock', 'virteng_params.yaml')
    pt_options.export_widgets_to_yaml('pretreatment_input', 'virteng_params.yaml', 'virteng_params.yaml')

    # Move into the pretreatment directory
    os.chdir('pretreatment_model/test/')
    
    # Run pretreatment code specifying location of input file
    path_to_input_file = '%s/virteng_params.yaml' % (parent_path)
    %run ptrun.py $path_to_input_file
    
    if pt_options.show_plots.value:
        %run postprocess.py 'out_\*.dat' exptdata_150C_1acid.dat

    os.chdir(parent_path)
    print('\nFinished Pretreatment')

    
def run_enzymatic_hydrolysis():
    print('\nRunning Enzymatic Hydrolysis Model')
    parent_path = os.getcwd()

    # Export the enzymatic hydrolysis options to a global yaml file
    eh_options.export_widgets_to_yaml('enzymatic_input', 'virteng_params.yaml', 'virteng_params.yaml')
    
    # Run the selected enzymatic hydrolysis model
    if eh_options.use_cfd.value:
        
        # Export the current state to a dictionary
        virteng_params = yaml_to_dict('virteng_params.yaml') 
        
        enzdata_replacements = {}
        enzdata_replacements['lmbde'] = virteng_params['enzymatic_input']['lambda_e']
        enzdata_replacements['fis0'] = virteng_params['enzymatic_input']['fis_0']
        enzdata_replacements['dt_ss'] = virteng_params['enzymatic_input']['t_final']
    
        os.chdir('EH_CFD/')

        write_file_with_replacements('enzdata', enzdata_replacements)
        
        if hpc_run:
            host_list = !srun hostname
            num_nodes = len(host_list)
            max_cores = int(36*num_nodes)
            
            !./run.sh $max_cores
        else:
            print('Cannot run EH_CFD without HPC resources.')
            print('$ ./run.sh $max_cores')
            print(os.getcwd())

    else:
        os.chdir('two_phase_batch_model/')
        path_to_input_file = '%s/virteng_params.yaml' % (parent_path)
        %run two_phase_batch_model.py $path_to_input_file
    
    os.chdir(parent_path)
    print('\nFinished Enzymatic Hydrolysis')

    
def run_bioreactor():
    print('\nRunning Bioreactor Model')
    parent_path = os.getcwd()

    # Export the bioreactor options to a global yaml file
    br_options.export_widgets_to_yaml('bioreactor_input', 'virteng_params.yaml', 'virteng_params.yaml')
    
    # Convert the current parameters file to a dictionary
    virteng_params = yaml_to_dict('virteng_params.yaml')        

    # Make changes to the fvOptions file based on replacement options
    fvOptions_replacements = {}
    for key, value in virteng_params['enzymatic_output'].items():
        fvOptions_replacements['double %s' % (key)] = value

    os.chdir('bioreactor/bubble_column/constant/')
    write_file_with_replacements('fvOptions', fvOptions_replacements)
    os.chdir(parent_path)
    
    # Make changes to the controlDict file based on replacement options
    controlDict_replacements = {}
    controlDict_replacements['endTime '] = virteng_params['bioreactor_input']['t_final']
    
    os.chdir('bioreactor/bubble_column/system/')
    write_file_with_replacements('controlDict', controlDict_replacements)
    os.chdir(parent_path)

    # Run the bioreactor model
    os.chdir('bioreactor/bubble_column/')
    if hpc_run:
        # call function to update ovOptions # fvOptions?
        !sbatch ofoamjob
    else:
        print('Cannot run bioreactor without HPC resources.')
        print('$ sbatch ofoamjob')
        print(os.getcwd())
        
    os.chdir(parent_path)
    print('\nFinished Bioreactor')


# Define a function to be executed each time the run button is pressed
def run_button_action(b):
    clear_output()
    display(run_button)
    
    # Run the pretreatment model
    run_pretreatment()
    
    # Run the enzymatic hydrolysis model
    run_enzymatic_hydrolysis()
    
    # Run the bioreactor model
    run_bioreactor()
    
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

```

```python

```
