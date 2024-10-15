import sys
import os
import subprocess
import numpy as np

from virteng.FileModifiers import write_file_with_replacements
from virteng.Utilities import check_dict_for_nans, dict_to_yaml, yaml_to_dict, print_dict
from virteng.WidgetFunctions import OptimizationWidget
from virteng.ModelsConnection import VE_params

root_path = os.path.abspath(os.path.join(os.path.dirname(__file__)))
cwd = os.getcwd()

def make_models_list(options_list, n_models=2, hpc_run=False):

    # Initialize models
    model1 = Model1(options_list[0], hpc_run)
    model2 = Model2(options_list[1])
    models_list = [model1, model2]
    return models_list

def make_output_names():
    return ['model1_out', 'model2_out']


class Model1:
    def __init__(self, model1_options, hpc_run=False):
        ''' 
        :param model1_options: (WidgetCollection) or (dict)
            A ``WidgetCollection`` object containing all of widgets used
            to solicit user input for feedstock properties 
            or dictionary with feedstock properties.
        '''
        self.ve= VE_params()
        self.ve.model1_in = {}
        if type(model1_options) is dict:
            self.freq_factor = model1_options['freq_factor']
            self.act_energy = model1_options['act_energy']
            self.temp = model1_options['temp']
        else:
            for widget_name, widget in model1_options.__dict__.items(): 
                if isinstance(widget, OptimizationWidget):
                    setattr(self, widget_name, widget.widget.value)
                else:
                    setattr(self, widget_name, widget.value)
        self.hpc_run = hpc_run
        self.model1_module_path = os.path.join(root_path, 'models', 'model1')
        sys.path.append(self.model1_module_path)

    ##############################################
    ### Properties
    ##############################################
    @property
    def freq_factor(self):
        return self.ve.model1_in['freq_factor']

    @freq_factor.setter
    def freq_factor(self, a):
        if not 0 <= a:
            raise ValueError(f"Value of var1 = {a} is negative")
        self.ve.model1_in['freq_factor'] = float(a)

    @property
    def act_energy(self):
        return self.ve.model1_in['act_energy']

    @act_energy.setter
    def act_energy(self, a):
        if not 0 <= a:
            raise ValueError(f"Value of var2 = {a} is negative")
        self.ve.model1_in['act_energy'] = float(a)

    @property
    def temp(self):
        return self.ve.model1_in['temp']

    @temp.setter
    def temp(self, a):
        if not 0 <= a:
            raise ValueError(f"Value of var3 = {a} is negative")
        self.ve.model1_in['temp'] = float(a)

    
    ##############################################
    ### run model
    ##############################################
    def run(self, verbose=False):

        if not self.hpc_run:
            from model1 import run_model1
            self.ve.model1_out = run_model1(self.ve) 
        else:
            os.chdir(self.model1_module_path)
            self.ve.write_to_file('ve_params.yml')
            command = 'python run_model1.py'
            subprocess.call(command.split(), text=True)
            self.ve = VE_params.load_from_file('ve_params.yml', verbose=False)
            os.chdir(root_path)
        
        if check_dict_for_nans(self.ve.model1_out):
            return True
        return False

class Model2:
    def __init__(self, model2_options):
        ''' 
        :param model1_options: (WidgetCollection) or (dict)
            A ``WidgetCollection`` object containing all of widgets used
            to solicit user input for feedstock properties 
            or dictionary with feedstock properties.
        '''
        self.ve= VE_params()
        self.ve.model2_in = {}
        if type(model2_options) is dict:
            self.time = model2_options['time']
        else:
            for widget_name, widget in model2_options.__dict__.items(): 
                if isinstance(widget, OptimizationWidget):
                    setattr(self, widget_name, widget.widget.value)
                else:
                    setattr(self, widget_name, widget.value)

        self.model2_module_path = os.path.join(root_path, 'models', 'model2')
        sys.path.append(self.model2_module_path)

    ##############################################
    ### Properties
    ##############################################
    @property
    def time(self):
        return self.ve.model2_in['time']

    @time.setter
    def time(self, a):
        if not 0 <= a <= 100:
            raise ValueError(f"Value {a} is outside allowed interval [0, 100]")
        self.ve.model2_in['time'] = float(a)


    
    ##############################################
    ### run model
    ##############################################
    def run(self, verbose=False):

        from model2 import run_model2
        self.ve.model2_out = run_model2(self.ve) 

        
        if check_dict_for_nans(self.ve.model1_out):
            return True
        return False