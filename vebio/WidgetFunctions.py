from ipywidgets import *

from vebio.Utilities import dict_to_yaml, yaml_to_dict

#================================================================

class WidgetCollection:
    def __init__(self):
        pass

    def display_all_widgets(self):
        
        # Set default viewing options
        description_width = 175
        widget_width = 350
        info_width = 300
        padding = 10
        
        # Define display options
        widget_style = {'description_width': '%dpx' % (description_width)}
        info_layout = {'margin': '0px 0px 0px %dpx' % (2*padding), 'width':'%dpx' % (info_width)}
        # box_layout = {'padding': '10px'}
        box_layout = {'padding': '0px %dpx 0px %dpx' % (padding, padding), 'align_items': 'center'}

        # For every widget
        for item in self.__dict__.items():

            widget_layout = {'width': '%dpx' % (widget_width)}

            # Extract the first object
            w = item[1]

            if hasattr(w, 'contains_sub_widgets'):
                w.lower.style = widget_style
                
                sub_width = int((widget_width - description_width)/2.0 - 2.0)
                
                w.lower.layout = {'width': '%dpx' % (description_width + sub_width)}
                w.upper.layout = {'width': '%dpx' % (sub_width)}
                
                myLabel = widgets.HTMLMath(
                    value = w.lower.description_tooltip,
                    layout = info_layout
                )
                
                w.box = HBox([w.lower, w.upper, myLabel], layout = box_layout)

            else:
                # Set this widget's style and layout
                w.style = widget_style

                if type(w) == Checkbox:
                    shift_amt = 0.88*(widget_width - description_width)
                    widget_layout.update({'padding': '0px 0px 0px %dpx ' % (shift_amt)})

                if hasattr(w, 'custom_layout'):
                    widget_layout.update(w.custom_layout)

                w.layout = widget_layout

                myLabel = widgets.HTMLMath(
                    value = w.description_tooltip,
                    layout = info_layout
                )

                # Organize this widget with more layout options
                w.box = HBox([w, myLabel], layout = box_layout)

            display(w.box)

    def export_widgets_to_dict(self, parent_name=None):
        
        #Start with a blank dictionary
        widget_dict = {}

        for widget_name, widget in self.__dict__.items():
            # Get the name and current state of each widget

            widget_value = widget.value

            if hasattr(widget, 'scaling_fn'):
                # print('pre-scaling value = %f' % (widget_value))
                widget_value = widget.scaling_fn(widget_value)
                # print('post-scaling value = %f' % (widget_value))

            # Create a dictionary with name : value pairs
            widget_dict['%s' % (widget_name)] = widget_value

        if parent_name is not None:
            widget_dict = {'%s' % (parent_name): widget_dict}

        return widget_dict

#================================================================

class ValueRangeWidget:
    def __init__(self, description, tooltip, bounds, init_vals, step_size):
        self.contains_sub_widgets = True
        
        self.lower = widgets.BoundedFloatText(
            value = init_vals[0],
            min = bounds[0],
            max = bounds[1],
            step = step_size,
            description = description,
            description_tooltip = tooltip,
            disabled = False
        )

        self.upper = widgets.BoundedFloatText(
            value = init_vals[1],
            min = bounds[0],
            max = bounds[1],
            step = step_size,
            disabled = False
        )

        def swap_range_values(change):
            lower_val_tmp = self.lower.value
            upper_val_tmp = self.upper.value

            if upper_val_tmp < lower_val_tmp:
                self.lower.value = upper_val_tmp
                self.upper.value = lower_val_tmp

        self.upper.observe(swap_range_values, names='value')
        self.lower.observe(swap_range_values, names='value')
