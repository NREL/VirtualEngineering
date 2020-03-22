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


```python
%reset -f
from ipywidgets import *
from IPython.display import HTML, clear_output
import yaml
import os
import time
from shutil import copyfile

#================================================================

# See if we're running on Eagle or on a laptop

hostname_output = !hostnamectl

hpc_test = [line for line in hostname_output if 'computer-server' in line]

if len(hpc_test) > 0:
    hpc_run = True
else:
    hpc_run = False
    print('It looks like you\'re running this notebook on a laptop.')
    print('Some features requiring HPC resources will be disabled.')

#================================================================
# I think it would be cleaner to have these function definitions in a separate file
# to be read in, JJS 3/22/20
def display_all_widgets(widget_collection):
    
    # Define display options
    widget_style = {'description_width': '200px'}
    widget_layout = {'width': '400px'}
    info_layout = {'margin': '0px 0px 0px 10px', 'width':'400px'}
    # box_layout = {'padding': '10px'}
    box_layout = {'padding': '10px', 'align_items': 'center'}

    # For every widget
    for item in widget_collection.__dict__.items():
        # Extract the first object
        w = item[1]

        # Set this widget's style and layout
        w.style = widget_style
        w.layout = widget_layout
        
        myLabel = widgets.HTMLMath(
            value = w.description_tooltip,
            layout = info_layout
        )
        
        # Organize this widget with more layout options
        w.box = HBox([w, myLabel], layout = box_layout)

        display(w.box)

#================================================================

def export_widgets_to_yaml(widget_collection, yaml_filename, merge_output_file=None):
    
        #Start with a blank dictionary
        widget_dict = {}

        for item in widget_collection.__dict__.items():
            # Get the name and current state of each widget
            widgetName = item[0]
            widgetValue = item[1].value
            
            # Create a dictionary with name : value pairs
            widget_dict['%s' % (widgetName)] = widgetValue
            
        if merge_output_file is not None:
            with open(merge_output_file) as fp:
                merge_dict = yaml.load(fp, Loader = yaml.FullLoader)
                
            widget_dict.update(merge_dict)
            
        # Dump the new dictionary into a yaml file
        with open(yaml_filename, 'w') as fp:
            yaml.dump(widget_dict, fp)
            
#================================================================

class blank_object:
    pass

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
fs_options = blank_object()

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
display_all_widgets(fs_options)

#================================================================

```

---

### 1. Pretreatment Operation

Set the options for the pretreatment operation below.

```python
#================================================================

# Create the collection of widgets
pt_options = blank_object()

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
    value = 2400,
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
display_all_widgets(pt_options)

#================================================================

```

---

### 2. Enzymatic Hydrolysis Operation

Set the options for the enzymatic hydrolysis operation using a two-phase reaction rate model below.


```python
#================================================================

# Create the collection of widgets
eh_options = blank_object()

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

# eh_options.dilution_strength = widgets.FloatSlider(
#     value = 0.5,
#     min = 0,
#     max = 1,
#     step = 0.1,
#     description = 'Dilution Strength',
#     continuous_update = False,
#     orientation = 'horizontal',
#     readout = True,
#     readout_format = '.1f',
#     description_tooltip = r'The efficacy of the dilution step: 0 means use the FIS$_0$ value from the end of the pretreatment step entirely, 1 means use the target FIS$_0$ value entirely'
# )


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
display_all_widgets(eh_options)

#================================================================

```

---

### 3. Bioreaction Operation

Set the options for the bubble column bioreaction operation below.


```python
#================================================================

# Create the collection of widgets
br_options = blank_object()

br_options.t_final = widgets.BoundedFloatText(
    value = 100.0,
    min = 1.0,
    max = 1e16,
    description = 'Final Time',
    description_tooltip = r'The total time of the simulation (h).  Must be $\geq 1$'
)

#================================================================

# Display the widgets
display_all_widgets(br_options)

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
    os.chdir('pretreatment_model/test/')
    
    # These print statements do not show for me until the run is complete. Something that should be fixed
    # at some point, JJS 3/22/20
    print('Running Pretreatment Model')
    export_widgets_to_yaml(fs_options, 'fs_input.yaml')
    export_widgets_to_yaml(pt_options, 'pt_input.yaml', 'fs_input.yaml')
    %run ptrun.py 'pt_input.yaml' 'pt_output.yaml'
    if pt_options.show_plots.value:
        %run postprocess.py 'out_\*.dat' exptdata_150C_1acid.dat
    copyfile('pt_output.yaml', '../../two_phase_batch_model/%s' % ('pt_to_eh_input.yaml'))
    print('\nFinished Pretreatment')
    
    # Run the enzymatic hydrolysis model
    os.chdir(parent_path)
    os.chdir('two_phase_batch_model/')
    
    print('\nRunning Enzymatic Hydrolysis Model')
    export_widgets_to_yaml(eh_options, 'eh_input.yaml', 'pt_to_eh_input.yaml')
    # output argument is not currently used, but it should be for passing info to the bioreaction 
    # unit operation, JJS 3/22/20
    %run two_phase_batch_model_fitting.py 'eh_input.yaml' 'eh_output.yaml'
    print('\nFinished Enzymatic Hydrolysis')
    
    # Run the bioreactor model
    os.chdir(parent_path)
    os.chdir('bioreactor/bubble_column/')
    
    print('\nRunning Bioreactor Model')
    if hpc_run:
        # call function to update ovOptions
        !sbatch ofoamjob
    else:
        print('Cannot run bioreactor without HPC resources.')
        print('$ sbatch ofoamjob')
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
