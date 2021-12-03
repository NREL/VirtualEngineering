from ipywidgets import *

from vebio.Utilities import dict_to_yaml, yaml_to_dict

#================================================================

class WidgetCollection:
    def __init__(self):
        pass

    def display_all_widgets(self):
        
        # Set default viewing options
        widget_width = 350
        description_width = 125
        html_width = 350
        padding = 5
        
        # Define display options
        default_widget_layout = {'width': '%dpx' % (widget_width)}
        widget_style = {'description_width': '%dpx' % (description_width)}
        html_layout = {'width':'%dpx' % (html_width), 'margin': '0px 0px 0px %dpx' % (2*padding)}
        box_layout = {'padding': '0px %dpx 0px %dpx' % (padding, padding), 'align_items': 'center'}

        # For every widget
        for widget_name, widget in self.__dict__.items():

            widget_layout = default_widget_layout.copy()

            if hasattr(widget, 'contains_sub_widgets'):
                widget.lower.style = widget_style
                
                sub_width = int((widget_width - description_width)/2.0 - 2.0)
                
                widget.lower.layout = {'width': '%dpx' % (description_width + sub_width)}
                widget.upper.layout = {'width': '%dpx' % (sub_width)}
                
                html_label = widgets.HTMLMath(
                    value = widget.lower.description_tooltip,
                    layout = html_layout
                )
                
                hbox = HBox([widget.lower, widget.upper, html_label], layout = box_layout)

            else:
                # Set this widget's style and layout
                widget.style = widget_style

                if type(widget) == Checkbox:
                    shift_amt = (widget_width - description_width) - 22
                    widget_layout.update({'padding': '0px 0px 0px %dpx ' % (shift_amt)})

                elif type(widget) == RadioButtons:
                    height = (len(widget.options)-2)*20 + 2*24
                    widget_layout.update({'height': '%dpx' % (height)})

                if hasattr(widget, 'custom_layout'):
                    widget_layout.update(widget.custom_layout)

                widget.layout = widget_layout

                html_label = widgets.HTMLMath(
                    value = widget.description_tooltip,
                    layout = html_layout
                )

                # Organize this widget with more layout options
                hbox = HBox([widget, html_label], layout = box_layout)

            display(hbox)

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
