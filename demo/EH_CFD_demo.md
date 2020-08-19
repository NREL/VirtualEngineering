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

```python
from ipywidgets import *
from IPython.display import HTML, clear_output
import os

import vebio.WidgetFunctions as wf
import vebio.FileModifiers as fm

eh_cfd_options = wf.WidgetCollection()

eh_cfd_options.lambda_e = widgets.BoundedFloatText(
    value = 0.03,
    max = 1.0,
    min = 0.0,
    step = 0.01,
    description = 'Enzymatic Load',
    description_tooltip = 'Ratio of the enzyme mass to the total\
    solution mass (kg/kg).  Must be in the range [0, 1]'
)
eh_cfd_options.lambda_e.alt_name = 'lmbde'

eh_cfd_options.fis_0_target = widgets.BoundedFloatText(
    value = 0.05,
    max = 1.0,
    min = 0.0,
    step = 0.01,
    description = r'FIS$_0$ Target',
    description_tooltip = 'The target value for initial fraction of \
    insoluble solids *after* dilution (kg/kg).  Must be in the range [0, 1]'
)
eh_cfd_options.fis_0_target.alt_name = 'fis0'

eh_cfd_options.t_final = widgets.BoundedFloatText(
    value = 560.0,
    min = 1.0,
    max = 1e16,
    description = 'Final Time',
    description_tooltip = r'The total time of the simulation (s).  Must be $\geq$ 1'
)
eh_cfd_options.t_final.alt_name = 'dt_ss'


# Display the widgets
eh_cfd_options.display_all_widgets()

```

```python
run_button = widgets.Button(
    description = 'Run EH CFD Model.',
    tooltip = 'Execute the CFD simulation of the enzymatic hydrolysis \
    model with the options specified above.',
    # margin: top, right, bottom, left
    layout =  {'width': '200px', 'margin': '25px 0px 100px 170px'}, 
    button_style = 'success'
)

#================================================================

# run_button_output = widgets.Output()
display(run_button)

def run_button_action(b):
    clear_output()
    display(run_button)
    
    # Export the current state to a dictionary
    widget_dict = eh_cfd_options.export_widgets_to_yaml()

    print(widget_dict)

    parent_path = os.getcwd()

    os.chdir('../EH_CFD/')

    # fm.write_file_with_replacements(filename, replacement_dictionary)
    fm.write_file_with_replacements('enzdata', widget_dict)

    os.chdir(parent_path)

run_button.on_click(run_button_action)

```

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
