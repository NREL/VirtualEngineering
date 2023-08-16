
# import os
# import sys
import scipy.optimize as opt
import numpy as np

# imports from vebio modules
from vebio.WidgetFunctions import OptimizationWidget
# from vebio.Utilities import  yaml_to_dict
from vebio.RunFunctions import Feedstock, Pretreatment, EnzymaticHydrolysis, Bioreactor, VE_params
# # add path for no-CFD EH model
# sys.path.append(os.path.join(notebookDir, "submodules/CEH_EmpiricalModel/"))

class Optimization:

    def __init__(self, fs_options, pt_options, eh_options, br_options, obj_widjet, hpc_run):

        self.hpc_run  = hpc_run
        
        self.ve = VE_params()
        
        # self.output_names = [key for key in self.ve.__dict__.keys() if "out" in key]
        self.output_names = ['pt_out', 'eh_out', 'br_out']
        self.output_name = obj_widjet.value[0]
        self.objective_name = obj_widjet.value[-1]
        self.n_models = self.define_n_models()

        # Initialize models
        self.FS_model = Feedstock(fs_options)
        self.PT_model = Pretreatment(pt_options, hpc_run)
        self.models_list = [self.FS_model, self.PT_model]
        if self.n_models > 1:
            assert eh_options.model_type.value != 'CFD Simulation'
            self.EH_model = EnzymaticHydrolysis(eh_options, hpc_run)
            self.models_list.append(self.EH_model)
        if self.n_models > 2:
            assert br_options.model_type.value == 'CFD Surrogate' # Do optimization only with surrogate
            self.BR_model = Bioreactor(br_options, hpc_run)
            self.models_list.append(self.BR_model)

        self.fn_evals = 0
        self.objective_scaling = 1.0
        
        # Find variable to be optimized and set initial values and bounds 
        self.x_0 = []
        self.var_names = []
        self.nice_var_names = []
        self.var_bounds = []
        self.var_real_bounds = []
        for wc in [fs_options, pt_options, eh_options]:
            for widget_name, widget in wc.__dict__.items():        
                if isinstance(widget, OptimizationWidget) and widget.is_control.value == True:
                    bounds = (widget.widget.min, widget.widget.max)
                    scaled_value = self.normalize(widget.widget.value, bounds)        
                    # current_val = widget.widget.value

                    # Here, we use the values of the controls scaled to the range [0, 1]
                    self.x_0.append(scaled_value)
                    self.var_names.append(widget_name)
                    self.nice_var_names.append(widget.widget.description)
                    self.var_real_bounds.append(bounds)
                    self.var_bounds.append((0.0, 1.0))

                    print('Optimizing %s.' % widget.widget.description)

        if len(self.x_0) == 0:
            raise ValueError('No controls have been specified, retry with >= 1 control variables.')

    @staticmethod
    def normalize(value, bounds):
        lb, ub = bounds
        return (value - lb)/(ub-lb)
    
    @staticmethod
    def scale_back(value, bounds):
        lb, ub = bounds
        return value*(ub-lb) + lb

    def scipy_minimize(self, objective_fn, method='SLSQP', opt_results_file='optimization_results.csv'):
        self.opt_results_file = opt_results_file
        # Write header for the outputfile
        with open(self.opt_results_file, 'w') as fp:
            fp.write('# Iteration, ')
            for control in self.nice_var_names:
                fp.write('%s, ' % (control))
            fp.write('Objective\n')
        # Minimization
        self.opt_result = opt.minimize(objective_fn, self.x_0, method=method,
                                       bounds=self.var_bounds, callback=self.opt_callback)
        return self.opt_result

    def define_n_models(self):
        if not self.output_name in self.output_names:
            raise ValueError(f"Error: Output dictionary '{self.output_name}' doesn't exist. Check the widget definition")
        n = self.output_names.index(self.output_name) + 1
        print(f'Objective "{self.objective_name}" is in {self.output_name}.')
        print(f'On each iteration running n={n} models\n')
        return n

    def run_models_with_new_values(self, dimensional_values, verbose=False):
        # Update the models with the latest values
        for var_name, value in zip(self.var_names, dimensional_values):
            for model in self.models_list:
                if hasattr(model, var_name):
                    setattr(model, var_name, value)
                    
        # Run models
        for model in self.models_list:
            flag_nan = model.run(verbose=verbose)
            if flag_nan:
                return np.nan
        # Read the outputs into a dictionary
        obj = getattr(self.ve, self.output_name)[self.objective_name]
        return obj

    def objective_function(self, free_variables):
        
        # Scale back to dimensional values
        dimensional_values = []
        for value, bounds in zip(free_variables, self.var_real_bounds):
            dimensional_values.append(self.scale_back(value, bounds))
        # print('dimensional_values:', dimensional_values)

        # We take the negative so the minimize function sees the correct orientation
        obj = -self.run_models_with_new_values(dimensional_values, verbose=False)        

        # Set objactive scaling to normalize objective function to -1 before iterations 
        if self.fn_evals == 0:
            print('\nBeginning Optimization')
            print('objective scaling:', self.objective_scaling)
            self.objective_scaling = -1.0/obj
            
        self.fn_evals += 1
        
        # Write iteration in the file
        with open(self.opt_results_file, 'a') as fp:
            fp.write('%d, ' % (self.fn_evals))
            for dv in dimensional_values:
                fp.write('%.15e, ' % (dv))
            fp.write('%.15e\n' % (obj))

        print('Iter = %3d: ' % (self.fn_evals), end='')
        for k, dv in enumerate(dimensional_values):
            print('%s = %12.9e, ' % (self.var_names[k], dv), end='')
        print('Objective = %12.9e' % (-obj))
        
        obj *= self.objective_scaling
        # print(f'Scaled objective: {obj}')
        
        return obj

    @staticmethod
    def opt_callback(free_variables):
        pass
        # print('Controls:', free_variables)


    def parameter_grid_sweep(self, nn, results_file='sweep_params.cvs', verbose=False):
        """ Sample input parameters on a grid and run simulations.

        :param nn: (int) The number of points to select across each value
        :param results_file: The filename to write sweep results including 
                             extension, defaults to 'sweep_params.cvs'
        """
        
        with open(results_file, 'w') as fp:
            str_names = ''
            for name in self.nice_var_names:
                str_names += name + ', '  
            fp.write(f'# Iteration, {str_names}{self.objective_name}\n')
        
        # Make parameter grid
        def uniform_grid(bounds, nn):
            C_tmp = np.linspace(bounds[0], bounds[1], nn + 1)
            return C_tmp[:-1] + (C_tmp[1] - C_tmp[0]) / 2

        grid_x = [uniform_grid(bounds, nn) for bounds in self.var_real_bounds]
        grid_mesh = np.meshgrid(*grid_x, indexing='ij')
        dimension = len(self.var_names)
        grid_ravel = np.empty((dimension, nn**dimension))
        for i in range(dimension):
            grid_ravel[i] = grid_mesh[i].ravel() 
        
        for i, dimensional_values in enumerate(grid_ravel.T):
            # print(f'Iteration #{i+1}, parameter values: {dimensional_values}')
            # Run models 
            obj = self.run_models_with_new_values(dimensional_values, verbose=verbose)
            # Write output
            with open(results_file, 'a') as fp:
                str_values = ''
                for value in dimensional_values:
                    str_values += f'{value:.9e}' + ', '  
                fp.write(f'{(i+1):.0f}, {str_values}{obj:.9e}\n')
        
        print(f'Finished {nn**2} forward runs.')
        
        print('\nFinished sweeps!')
