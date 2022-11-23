from ipywidgets import *

from vebio.Utilities import dict_to_yaml, yaml_to_dict

#================================================================

class WidgetCollection:
    """A ``WidgetCollection`` object collects any number of different iPyWidgets.

    This object can be used to organize all the widgets pertaining to a single
    unit model or step of the overall conversion process.  New widgets should be
    accumulated after initialization using an approach like::


        fs_options = wf.WidgetCollection()
        fs_options.xylan_solid_fraction = widgets.BoundedFloatText(...)
        fs_options.glucan_solid_fraction = widgets.BoundedFloatText(...)
        fs_options.initial_porosity = widgets.Checkbox(...)

    where widget options and values can be accessed using, for example,
    ``fs_options.initial_porosity.value``

    """
    def __init__(self):

        pass

    def display_all_widgets(self):
        """Displays all the collected widgets in the Notebook.

        This method displays all the widgets collected in the object
        to the Jupyter Notebook interface.

        Args:
            None

        Returns:
            None

        """

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

            if hasattr(widget, 'custom_display'):
                widget.custom_display()

            else:
                # Set this widget's style and layout
                widget.style = widget_style

                if type(widget) == Checkbox:
                    shift_amt = (widget_width - description_width) - 22
                    widget_layout.update({'padding': '0px 0px 0px %dpx ' % (shift_amt)})

                elif type(widget) == RadioButtons:
                    height = (len(widget.options)-2)*20 + 2*24
                    widget_layout.update({'height': '%dpx' % (height)})

                    desc_widg_space = 5
                    desc_width = description_width + desc_widg_space
                    widget_layout.update({'width': f'{widget_width-desc_width:.0f}px'})

                    widget.desc_label = widgets.HTMLMath(
                        value = widget.description,
                        layout = {'width': f'{desc_width:.0f}px',
                        'display': 'flex',
                        'justify_content': 'flex-end',
                        'padding': '0px 5px 0px 0px'}
                    )

                    widget.description = ""

                if hasattr(widget, 'custom_layout'):
                    widget_layout.update(widget.custom_layout)

                widget.layout = widget_layout

                html_label = widgets.HTMLMath(
                    value = widget.description_tooltip,
                    layout = html_layout
                )

                # Organize this widget with more layout options
                if hasattr(widget, 'desc_label'):
                    hbox = HBox([widget.desc_label, widget, html_label], layout = box_layout)
                else:
                    hbox = HBox([widget, html_label], layout = box_layout)

                display(hbox)

    def export_widgets_to_dict(self, parent_name=None):
        """Store all widget values in dictionary.

        This method allows the values of each widget to be saved
        in a dictionary with the pattern ``{"widget_1_name": widget_1_value, ...}``.
        If the widget was created with a scaling function, the scaled
        version of the value is calculated and stored in the dictionary.
        This scaled value is the one that will be referenced by subsequent operations
        accessing the VE parameter file.

        Args:
            parent_name (str, optional):
                At the end of the dictionary creation, all the ``name: value``
                entries can be nested under a single "parent" keyword.  For example::

                    {parent_name: {"widget_1_name": widget_1_value,
                                   "widget_2_name": widget_2_value}}

                Defaults to ``None``, i.e., do not nest under a parent name.

        Returns:
            dict:
                The name and value of each widget
                collected in a Python dictionary.
            
        """

        #Start with a blank dictionary
        widget_dict = {}

        for widget_name, widget in self.__dict__.items():
            # Get the name and current state of each widget

            if isinstance(widget, OptimizationWidget):
                widget_value = widget.widget.value

            else:
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
    """A pair of linked ``BoundedFloatText`` widgets for value ranges.

    This is a custom widget composed of two linked ``BoundedFloatText`` widgets
    intended to enable the easy solicitation of a range of values from the user.

        
    """

    def __init__(self, description, tooltip, bounds, init_vals, step_size):

        """Initialize the linked input fields.

        These linked widgets will be displayed side-by-side for intuitive
        entry of [lower_bound, upper_bound] information, the behavior of the
        individual fields and their labeling/description is provided by the user.

        Args:
            description (string):
                A short name or label for the the widget, e.g., `Porosity Range`
            tooltip (string):
                A longer, potentially multi-line explanation describing the 
                values that are intended to be entered, background information,
                units, suggested values, etc., should all be communicated here.
            bounds (list(float)):
                A two-element list where ``bounds[0]`` is the lower bound to enforce
                on both fields and ``bounds[1]`` is the upper bound to enforce on
                both fields.
            init_vals (list(float)):
                A two-element list where ``init_vals[0]`` is the initial value of the
                lower bound and ``init_vals[1]`` is the initial value of the upper bound.
            step_size (float):
                The change in value produced by clicking the increment or decrement
                buttons in the ``BoundedFloatText`` field.

        Returns:
            None

        """

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

    def custom_display(self):
        """Displays a value range widget.

        This displays a value range widget including the formatting
        of the two fields side by side with a single label.

        Args:
            None
            
        Returns:
            None

        """

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

        self.lower.style = widget_style
        
        sub_width = int((widget_width - description_width)/2.0 - 2.0)
        
        self.lower.layout = {'width': '%dpx' % (description_width + sub_width)}
        self.upper.layout = {'width': '%dpx' % (sub_width)}
        
        self.html_label = widgets.HTMLMath(
            value = self.lower.description_tooltip,
            layout = html_layout
        )
        
        hbox = HBox([self.lower, self.upper, self.html_label], layout = box_layout)

        display(hbox)



#================================================================

class OptimizationWidget:
    """A wrapper around ipywidgets with options for optimization problems.

    This custom widget wraps around a specified ipywidget and includes
    checkboxes to specify whether or not this widget should be used as
    a control value or an objective value.

    """

    def __init__(self, widget_name, widget_options):
        
        """Initialize an optimization widget.

        This creates an optimization widget of type ``widget_name``
        which must be a valid ipywidget type. ``widget_options`` determines
        the set up and properties of this widget.

        Args:
            widget_name (string):
                The name of the ipywidget to create, e.g., `BoundedFloatText`
            widget_options (dict):
                A dictionary containing all the setup values for a widget of type
                ``widget_name``. For example::

                    widget_options = {'value': 0.25,
                                      'min': -1.0,
                                      'max': 1.0,
                                      'description': Short Description
                                      'description_tooltip': A long description...}

        Returns:
            None

        """

        self.widget = getattr(widgets, widget_name)(**widget_options)

        self.is_control = widgets.Checkbox(value = False,
                                           description = 'Control',
                                           disabled = False)

        def control_button(change):
            return

        self.is_control.observe(control_button, names='value')


    def custom_display(self):
        """Displays an optimization widget.

        This displays an optimization widget including the formatting
        of the boxes used to toggle control/objective specification.

        Args:
            None

        Returns:
            None

        """

        widget_width = 350
        description_width = 125

        checkbox_width = 100
        checkbox_description_width = 0

        html_width = 350
        padding = 5
        
        # Define display options
        widget_layout = {'width': f'{widget_width:.0f}px'}
        widget_style = {'description_width': f'{description_width:.0f}px'}

        checkbox_layout = {'width': f'{checkbox_width:.0f}px'}
        checkbox_style = {'description_width': f'{checkbox_description_width:.0f}px'}

        html_layout = {'width': f'{html_width:.0f}px',
                       'margin': f'0px 0px 0px {2.0*padding:.0f}px'}

        box_layout = {'padding': f'0px {padding:.0f}px 0px {padding:.0f}px',
                      'align_items': 'center'}

        self.widget.layout = widget_layout
        self.widget.style = widget_style

        self.is_control.layout = checkbox_layout
        self.is_control.style = checkbox_style
        
        self.html_label = widgets.HTMLMath(
            value = self.widget.description_tooltip,
            layout = html_layout
        )
        
        hbox = HBox([self.widget, self.is_control,  self.html_label], layout = box_layout)


        display(hbox)
