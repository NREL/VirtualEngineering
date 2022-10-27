
import os
import sys
import scipy.optimize as opt

# imports from vebio modules
from vebio.WidgetFunctions import OptimizationWidget
from vebio.Utilities import  yaml_to_dict
from vebio.RunFunctions import Feedstock, Pretreatment, EnzymaticHydrolysis, Bioreactor
# # add path for no-CFD EH model
# sys.path.append(os.path.join(notebookDir, "submodules/CEH_EmpiricalModel/"))

class Optimization:

    def __init__(self, fs_options, pt_options, eh_options, br_options, hpc_run, notebookDir, 
                params_filename='virteng_params_optimization.yaml', 
                opt_results_file='optimization_results.csv'):

        self.hpc_run  = hpc_run
        self.notebookDir = notebookDir
        self.params_filename = params_filename
        self.opt_results_file = opt_results_file 

        self.FS_model = Feedstock(params_filename, fs_options)
        self.PT_model = Pretreatment(notebookDir, params_filename, pt_options)
        self.EH_model = EnzymaticHydrolysis(notebookDir, params_filename, eh_options, hpc_run)
        self.BR_model = Bioreactor(notebookDir, params_filename, br_options, hpc_run)

        # Do optimization only with surrogates
        assert br_options.model_type.value == 'CFD Surrogate'
        assert eh_options.model_type.value == 'CFD Surrogate'

        self.fn_evals = 0
        self.objective_scaling = 1.0
        
        # Find variable to be optimized and set initial values and bounds 
        self.x_0 = []
        self.var_names = []
        self.var_bounds = []
        self.var_real_bounds = []
        for wc in [fs_options, pt_options, eh_options]:
            for widget_name, widget in wc.__dict__.items():        
                if isinstance(widget, OptimizationWidget) and widget.is_control.value == True:
                    print('Optimizing %s.' % widget_name)
                    
                    bounds = (widget.widget.min, widget.widget.max)
                    scaled_value = self.normalize(widget.widget.value, bounds)        
                    # current_val = widget.widget.value

                    # Here, we use the values of the controls scaled to the range [0, 1]
                    self.x_0.append(scaled_value)
                    self.var_names.append(widget_name)
                    self.var_real_bounds.append(bounds)
                    self.var_bounds.append((0.0, 1.0))

        if len(self.x_0) == 0:
            raise ValueError('No controls have been specified, retry with >= 1 control variables.')

        # Write header for the outputfile
        with open(self.opt_results_file, 'w') as fp:
            fp.write('# Iteration, ')
            for control in self.var_names:
                fp.write('%s, ' % (control))
            fp.write('Objective\n')

    @staticmethod
    def normalize(value, bounds):
        lb, ub = bounds
        return (value - lb)/(ub-lb)
    
    @staticmethod
    def scale_back(value, bounds):
        lb, ub = bounds
        return value*(ub-lb) + lb

    def scipy_minimize(self, objective_fn, method='SLSQP'):
        self.opt_result = opt.minimize(objective_fn, self.x_0,
                                  method=method,
                                  bounds=self.var_bounds,
                                  callback=self.opt_callback)
        return self.opt_result

    def objective_function(self, free_variables):
        
        # Scale back to dimensional values
        dimensional_values = []
        for value, bounds in zip(free_variables, self.var_real_bounds):
            dimensional_values.append(self.scale_back(value, bounds))

        # Update the models with the latest values
        if self.fn_evals != 0:
            for var_name, value in zip(self.var_names, dimensional_values):
                for model in [self.FS_model, self.PT_model, self.EH_model]:
                    if hasattr(model, var_name):
                        model.update_values(**{var_name: value})
                    
        # Set global paths and files for communication between operations
        os.chdir(self.notebookDir)
        
        # Turn off printed outputs from unit operations
        v_flag = (self.fn_evals == 0)
        
        self.PT_model.run(verbose=v_flag)          # Running the pretreatment model
        self.EH_model.run(verbose=v_flag)          # Run the enzymatic hydrolysis model
        self.BR_model.run(verbose=v_flag)          # Run the bioreactor model
        
        # Read the outputs into a dictionary
        output_dict = yaml_to_dict(self.params_filename)
        
        # The objective function in this case is OUR, we take the negative
        # so the minimize function sees the correct orientation
        obj = -output_dict['bioreactor_output']['our']        

        if self.fn_evals == 0:
            print('Beginning Optimization')
            self.objective_scaling = -1.0/obj
            
        self.fn_evals += 1
        
        # Write iteration in the file
        with open(self.opt_results_file, 'a') as fp:
            fp.write('%d, ' % (self.fn_evals))
            for dv in dimensional_values:
                fp.write('%.15e, ' % (dv))
            fp.write('%.15e\n' % (obj))

        print('\nIter = %3d: ' % (self.fn_evals), end='')
        for k, dv in enumerate(dimensional_values):
            print('%s = %12.9e, ' % (self.var_names[k], dv), end='')
        print('Objective = %12.9e' % (obj))
        
        obj *= self.objective_scaling
        print(f'Objactive scaled: {obj}')
        
        return obj

    @staticmethod
    def opt_callback(free_variables):
        pass
        # print('Controls:', free_variables)