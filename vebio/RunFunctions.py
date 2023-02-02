import sys
import os
import contextlib
import subprocess
import glob

import numpy as np


from vebio.FileModifiers import write_file_with_replacements
from vebio.Utilities import yaml_to_dict, dict_to_yaml, check_dict_for_nans

root_path = os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir))


class VE_params(object):
    ''' This  class is used for storing Virtual Engineering parameters 
        so they can be accesed from any model. It uses the Borg pattern. 
        The Borg pattern (also known as the Monostate pattern) is a way to
        implement singleton behavior, but instead of having only one instance
        of a class, there are multiple instances that share the same state. In
        other words, the focus is on sharing state instead of sharing instance.
    '''

    __shared_state = {}

    def __init__(self):
        self.__dict__ = self.__shared_state

    def __str__(self):
        return self.__shared_state

    # TODO:
    # def from_file(params_filename):

    # def write_to_file():


class Feedstock:
    def __init__(self, fs_options):
        """Through the ``fs_options`` widgets, the user controls the following
            values:

                * The initial fraction of solids due to xylan (X_X)
                * The initial fraction of solids due to glucan (X_G)
                * The initial porous fraction of the biomass particles

        :param fs_options: (WidgetCollection)
            A ``WidgetCollection`` object containing all of widgets used
            to solicit user input for feedstock properties.
        """
        self.ve= VE_params()
        self.ve.feedstock = {}
        self.xylan_solid_fraction = fs_options.xylan_solid_fraction.value
        self.glucan_solid_fraction = fs_options.glucan_solid_fraction.value
        self.initial_porosity = fs_options.initial_porosity.widget.value

    ##############################################
    ### Properties
    ##############################################
    @property
    def xylan_solid_fraction(self):
        return self.ve.feedstock['xylan_solid_fraction']

    @xylan_solid_fraction.setter
    def xylan_solid_fraction(self, a):
        if not 0 <= a <= 1:
            raise ValueError(f"Value {a} is outside allowed interval [0, 1]")
        self.ve.feedstock['xylan_solid_fraction'] = float(a)

    @property
    def glucan_solid_fraction(self):
        return self.ve.feedstock['glucan_solid_fraction']

    @glucan_solid_fraction.setter
    def glucan_solid_fraction(self, a):
        if not 0 <= a <= 1:
            raise ValueError(f"Value {a} is outside allowed interval [0, 1]")
        self.ve.feedstock['glucan_solid_fraction'] = float(a)

    @property
    def initial_porosity(self):
        return self.ve.feedstock['initial_porosity']

    @initial_porosity.setter
    def initial_porosity(self, a):
        if not 0 < a < 1:
            raise ValueError(f"Value {a} is outside allowed interval (0, 1)")
        self.ve.feedstock['initial_porosity'] = float(a)

    def run(self):
        return False


###################################################################################
####
####        PRETREATMENT
####
##################################################################################
class Pretreatment:

    def __init__(self, pt_options):
        """ Through the ``pt_options`` widgets, the user controls the following
            values:

                * Acid Loading (float)
                * Steam Temperature (float)
                * Initial FIS_0 (float)
                * Final Time (float)
                * Show plots (bool)

        :param pt_options: (WidgetCollection) or (dict)
            A ``WidgetCollection`` object containing all of widgets used
            to solicit user input for pretreatment properties.
            or 
            A dictionary containing all input values
        """

        print('Initializing Pretreatment Model')

        self.show_plots = pt_options.show_plots.value 

        self.ve = VE_params()
        self.ve.pt_in = {}
        self.initial_acid_conc = pt_options.initial_acid_conc.widget.value
        self.steam_temperature = pt_options.steam_temperature.widget.value
        self.initial_solid_fraction = pt_options.initial_solid_fraction.widget.value
        self.final_time = pt_options.final_time.widget.value

        pt_module_path = os.path.join(root_path,'pretreatment_model')
        sys.path.append(os.path.join(pt_module_path, 'fenics'))

    ##############################################
    ### Properties
    ##############################################
    @property
    def initial_acid_conc(self):
        return self.ve.pt_in['initial_acid_conc']

    @initial_acid_conc.setter
    def initial_acid_conc(self, a):
        if not 0 <= a <= 1:
            raise ValueError(f"Value {a} is outside allowed interval [0.0, 1.0]")
        self.ve.pt_in['initial_acid_conc'] = float(a)

    @property
    def steam_temperature(self):
        return self.ve.pt_in['steam_temperature'] - 273.15

    @steam_temperature.setter
    def steam_temperature(self, a):
        if not 3.8 <= a <= 250.3:
            raise ValueError(f"Value {a} is outside allowed interval [3.8, 250.3]")
        self.ve.pt_in['steam_temperature'] = float(a) + 273.15 # Conversion from celsius to kelvin

    @property
    def initial_solid_fraction(self):
        return self.ve.pt_in['initial_solid_fraction']

    @initial_solid_fraction.setter
    def initial_solid_fraction(self, a):
        if not 0 < a < 1:
            raise ValueError(f"Value {a} is outside allowed interval (0, 1)")
        self.ve.pt_in['initial_solid_fraction'] = float(a)

    @property
    def final_time(self):
        return self.ve.pt_in['final_time'] / 60

    @final_time.setter
    def final_time(self, a):
        if not 1 <= a <= 1440:
            raise ValueError(f"Value {a} is outside allowed interval [1, 1440]")
        self.ve.pt_in['final_time'] = 60 * float(a)

    ##############################################
    #
    ##############################################
    def run(self, verbose=True, show_plots=None):
        """Run pretreatment code specified in pretreatment_model/test/ptrun.py

        :param verbose: (bool, optional) 
            Option to show print messages from executed file, default True.
        :param verbose: (bool, optional) 
            Option to show plots, default True.
        """
        print('\nRunning Pretreatment')
    
        if show_plots is None:
            show_plots = self.show_plots
        from run_pretreatment import run_pt
        self.ve.pt_out = run_pt(self.ve, verbose, show_plots)

        print('Finished Pretreatment')
        if check_dict_for_nans(self.ve.pt_out):
            return True
        return False

###################################################################################
####
####        ENZYMATIC HYDROLYSIS
####
##################################################################################
class EnzymaticHydrolysis:
    def __init__(self, eh_options, hpc_run):
        """ Initialize enzymatic hydrolysis class. Three 
            distinct variants are included in the virtual engineering code:
            (1) a two-phase model which makes a well-mixed assumption, (2)
            a pre-trained surrogate model informed from CFD runs, and (3) 
            the CFD simulation itself, where option (3) is accessible only
            with ``hpc_run=True``. The default unit operation is the 
            surrogate model.

            Through the ``eh_options`` widgets, the user controls the following
            values:

                * Model Type
                * Enzymatic Load (float)
                * FIS_0 Target (float)
                * Final Time (float)
                * Show plots (bool)

        :param notebookDir: (str):
            The path to the Jupyter Notebook, used to specify the location
            of the input file and reset the working directory after this operation
            is finished.
        :param eh_options: (WidgetCollection):
            A ``WidgetCollection`` object containing all of widgets used
            to solicit user input for enzymatic hydrolysis properties.
        :param hpc_run: (bool)
            A flag indicating whether or not the Notebook is being
            run on HPC resources, enable CFD only if True.
        """

        print('Initializing Enzymatic Hydrolysis Model')

        self.hpc_run = hpc_run
        self.show_plots = eh_options.show_plots.value

        self.ve = VE_params()
        self.ve.eh_in = {}
        # EH input parameters
        self.lambda_e = eh_options.lambda_e.widget.value  # Conversion from mg/g to kg/kg
        self.fis_0 = eh_options.fis_0.value
        self.t_final = eh_options.t_final.value
        self.model_type = eh_options.model_type.value # running select_run_function() inside
        
    ##############################################
    ### Properties
    ##############################################
    @property
    def lambda_e(self):
        return self.ve.eh_in['lambda_e'] * 1000 # Conversion from kg/kg to mg/g 

    @lambda_e.setter
    def lambda_e(self, a):
        if not 0 <= a <= 1000:
            raise ValueError(f"Value {a} is outside allowed interval [0, 1000]")
        self.ve.eh_in['lambda_e'] = float(a) / 1000 # Conversion from mg/g to kg/kg
        # self.input2yaml(rewrite=True)

    @property
    def fis_0(self):
        return self.ve.eh_in['fis_0']

    @fis_0.setter
    def fis_0(self, a):
        if not 0 <= a <= 1:
            raise ValueError(f"Value {a} is outside allowed interval [0, 1]")
        self.ve.eh_in['fis_0'] = float(a)
        # self.input2yaml(rewrite=True)

    @property
    def t_final(self):
        return self.ve.eh_in['t_final']

    @t_final.setter
    def t_final(self, a):
        if not 1 <= a <= 24:
            raise ValueError(f"Value {a} is outside allowed interval [1, 24]")
        self.ve.eh_in['t_final'] = float(a)

    @property
    def model_type(self):
        return self.ve.eh_in['model_type']

    @model_type.setter
    def model_type(self, a):
        if not a in ['CFD Simulation', "CFD Surrogate", 'Lignocellulose Model']:
            raise ValueError("Invalid value. Allowed options: 'CFD Simulation', 'CFD Surrogate', 'Lignocellulose Model'")
        self.ve.eh_in['model_type']= a
        self.select_run_function()
    ##############################################
    #
    ##############################################

    def select_run_function(self):
        # selected enzymatic hydrolysis model
        if self.model_type == 'CFD Simulation':
            assert self.hpc_run, f'Cannot run EH_CFD without HPC resources. \n {os.getcwd()}'
            self.run = self.run_eh_cfd_simulation
            # TODO: add import path 
        elif self.model_type == "CFD Surrogate":
            self.run = self.run_eh_cfd_surrogate
            eh_module_path = os.path.join(root_path,'EH_OpenFOAM', 'EH_surrogate')
            if not eh_module_path in sys.path:
                sys.path.append(eh_module_path)
                
        elif self.model_type == 'Lignocellulose Model':
            self.run = self.run_eh_lignocellulose_model
            eh_module_path = os.path.join(root_path ,'two_phase_batch_model')
            if not eh_module_path in sys.path:
                sys.path.append(eh_module_path)

    def get_globalVars(self):
        """ Prepare input values for EH CFD operation

        :return: globalVar dictionary
        """
        globalVars = {}
        globalVars['fis0'] = self.fix_0
        globalVars['xG0'] = self.ve.pt_out['X_G']
        globalVars['xX0'] = self.ve.pt_out['X_X']
        globalVars['XL0'] = 1.0 - globalVars['xG0'] - globalVars['xX0']
        globalVars['yF0'] = 0.2 + 0.6*self.ve.pt_out['conv']
        globalVars['lmbdE'] = self.lambda_e
        globalVars['rhog0'] = 0.0
        self.dilution_factor = self.fis_0/self.ve.pt_out['fis_0']
        globalVars['rhox0'] = self.ve.pt_out['rho_x'] * self.dilution_factor
        globalVars['rhosl0'] = 0.0

        return globalVars

    def run_eh_cfd_simulation(self, verbose=True):
        
        print('\nRunning Enzymatic Hydrolysis Model')
        globalVars = self.get_globalVars()
        write_file_with_replacements('constant/globalVars', globalVars)

        # Get reaction_update_time, fluid_update_time, and fluid_steadystate_time
        # in order to convert the user-specified t_final into the endTime definition
        # expected by the OpenFOAM simulation
        reaction_update_time = 1.0
        fluid_update_time = 250.0
        fluid_steadystate_time = 400.0

        with open('constant/EHProperties', 'r') as fp:
            for line in fp:
                if '#' not in line:
                    if 'reaction_update_time' in line:
                        reaction_update_time = float(line.split(']')[-1].split(';')[0])
                    elif 'fluid_update_time' in line:
                        fluid_update_time = float(line.split(']')[-1].split(';')[0])
                    elif 'fluid_steadystate_time' in line:
                        fluid_steadystate_time = float(line.split(']')[-1].split(';')[0])

        controlDict = {}
        fintime = fluid_steadystate_time + (self.t_final/reaction_update_time + 1.0)*fluid_update_time
        controlDict['endTime'] = fintime

        write_file_with_replacements('system/controlDict', controlDict)

        '''
        import numpy as np
        import subprocess
        import os
        
        def check_queue(username, jobname):
            command = 'squeue -u %s -t R,PD -n %s' % (username, jobname)
            out = subprocess.run(command.split(), capture_output=True, text=True)
            print(out.stdout)
            
            if username in out.stdout:
                # Job is running, do nothing
                print('Job is already running')
                job_id = out.stdout.strip().split('\\n')[-1].split()[0]
                
            else:
                # Job is not running, submit it
                command = 'sbatch --job-name=%s dummy_job.sbatch' % (jobname)
                out = subprocess.run(command.split(), capture_output=True, text=True)
                print(out.stdout)
                job_id = out.stdout.strip().split()[-1]
                
                with open('job_history.csv', 'a') as fp:
                    fp.write('%s\\n' % (job_id))
                
            print(job_id, len(job_id))
                
        username = os.environ['USER']
        
        check_queue(username, 'dummy_job')
        '''
        
        # command = "srun hostname"
        # host_list = subprocess.run(command.split(), capture_output=True).stdout.decode()
        # num_nodes = len(host_list)
        # max_cores = int(36*num_nodes)
        
        # Fill output dict with nans, so Bioreactor know we are still running
        self.ve.eh_out = dict(zip(['rho_g', 'rho_x', 'rho_sl', 'rho_f'], [np.nan]*4))

        username = os.environ['USER']
        jobname = 'eh_cfd'

        # Check the queue
        command = 'squeue -u %s -t R,PD -n %s' % (username, jobname)
        out = subprocess.run(command.split(), capture_output=True, text=True)

        if username in out.stdout:
            # Job is running, do nothing
            print('EH CFD job is already queued.')
            print(out.stdout)
            job_id = out.stdout.strip().split('\n')[-1].split()[0]
            # TODO: there is no use_previous_output widget in notebook, comment for now 
        #     if eh_options.use_previous_output.value:
        #         print('Using outputs from most recent finished simulation.')
        #         integrated_quantities = np.genfromtxt('old_integrated_quantities.dat') # mol/L
        #         output_dict = {}
        #         output_dict['rho_g'] = float(integrated_quantities[-1, -3])
        #         output_dict['rho_x'] = float(integrated_quantities[-1, -2])
        #         output_dict['rho_sl'] = float(integrated_quantities[-1, -1])
        #         output_dict['rho_f'] = float(ve.pt_out['rho_f']*self.dilution_factor)
        #         print('Success.')
        else:
            # Job is not running, submit it
            print('Submitting EH CFD job.')
            command = 'sbatch --job-name=%s ofoamjob' % (jobname)
            out = subprocess.run(command.split(), capture_output=True, text=True)
            print(out.stdout)
            job_id = out.stdout.strip().split()[-1]

            with open('job_history.csv', 'a') as fp:
                fp.write('%s\n' % (job_id))

        print('Job ID = %s' % (job_id))
       
        # Prepare output values from EH CFD operations
        # FIXME: rho_g should be value taken from CFD output - should be good now, please check, JJS 9/15/20
        
        # per Hari's email, glucose concentration in fort.44 is mol/L
        # integrated_quantities = np.genfromtxt('integrated_quantities.dat', delimiter=' ') # mol/L

        '''
        This code represents the conversion that used to be necessary for the NEK 5000 simulation
        outputs, it's preserved here for reference but shouldn't be necessary for the new
        OpenFOAM version of EH.  Although it still may be necessary to calculate a version of
        dilution_factor_final and use it to scale the final four output values.


        rho_g_final = float(c_g_output[-1, 1])*180 # g/L
        
        # back-calculate fis from conversion value
        conv_output = np.genfromtxt('fort.42')
        conversion = float(conv_output[-1,1])
        fis_final = ve_params['enzymatic_input']['fis_0']*(1 - conversion)
        ## if have non-glucan solids, e.g. lignin, then the calculation will be:
        # fis = fis_0*(1 - XG0*conversion)
        ## where XG0 is initial fraction of solids that is glucan
        
        # this dilution calculation is not correct and needs fixing, JJS 9/15/20
        #dilution_factor_final = fis_final/ve_params['enzymatic_input']['fis_0']
        
        # FIXME: dilution_factor_final should be (fis_final)/(ve_params['enzymatic_input']['fis_0'])
        # where fis_final is taken from CFD output
        dilution_factor_final = 1.0
        rho_x_final = rho_x0*dilution_factor_final
        rho_f_final = rho_f0*dilution_factor_final
        '''
        
        os.chdir(self.notebookDir)
        print('Finished Enzymatic Hydrolysis')
        if check_dict_for_nans(self.ve.eh_out):
            return True
        return False

    def run_eh_cfd_surrogate(self, verbose=True):

        print('\nRunning Enzymatic Hydrolysis Model')
        
        from EH_surrogate import run_eh
        self.ve.eh_out = run_eh(self.ve)        
        
        print('Finished Enzymatic Hydrolysis')
        if check_dict_for_nans(self.ve.eh_out):
            return True
        return False

    def run_eh_lignocellulose_model(self, verbose=True):
        
        print('\nRunning Enzymatic Hydrolysis Model')
        # Commenting out cellulose-only two-phase model to use lignocellulose
        # model, just in case we want to switch back or make both an
        # option. The lignocellulose model is superior.
        #run_script("two_phase_batch_model.py", path_to_input_file, verbose=verbose)
        
        from driver_batch_lignocell_EH_VE import run_eh_lingocell
        self.ve.eh_out = run_eh_lingocell(self.ve, self.show_plots)
        print('Finished Enzymatic Hydrolysis')
        if check_dict_for_nans(self.ve.eh_out):
            return True
        return False

###################################################################################
####
####        BIOREACTOR
####
##################################################################################
class Bioreactor:
    def __init__(self, br_options, hpc_run):
        """ Initialize the aerobic bioreaction operation using 
            user-specified properties. Two distinct models exist: (1) a 
            pre-trained surrogate model informed from CFD runs and (2) the
            full CFD simulation itself where option (2) is accessible only
            with ``hpc_run=True``.  The default option is the surrogate model.

            Through the ``br_options`` widgets, the user controls the following
            values:

                * Model Type
                * Final Time (float)

        :param notebookDir: (str):
            The path to the Jupyter Notebook, used to specify the location
            of the input file and reset the working directory after this operation
            is finished.
        :param br_options: (WidgetCollection):
            A ``WidgetCollection`` object containing all of widgets used
            to solicit user input for bioreaction properties.
        :param hpc_run: (bool):
            A flag indicating whether or not the Notebook is being
            run on HPC resources, enable CFD only if True.
        """

        print('Initializing Bioreactor Model')
        self.hpc_run = hpc_run

        self.ve = VE_params()
        self.ve.br_in = {}

        self.br_module_path = os.path.join(root_path,'bioreactor', 'bubble_column')
        # Bioreactor input parameters
        self.model_type = br_options.model_type.value # running select_run_function() inside
        self.gas_velocity = br_options.gas_velocity.value
        self.column_height = br_options.column_height.value
        self.column_diameter = br_options.column_diameter.value
        self.bubble_diameter = br_options.bubble_diameter.value
        self.t_final = br_options.t_final.value

    ##############################################
    ### Properties
    ##############################################
    @property
    def gas_velocity(self):
        return self.ve.br_in['gas_velocity']

    @gas_velocity.setter
    def gas_velocity(self, a):
        if not 0.01 <= a <=0.1:
            raise ValueError(f"Value {a} is outside allowed interval [1, 1e16]")
        self.ve.br_in['gas_velocity'] = float(a)

    @property
    def column_height(self):
        return self.ve.br_in['column_height']

    @column_height.setter
    def column_height(self, a):
        if not 10 <= a <= 50:
            raise ValueError(f"Value {a} is outside allowed interval [1, 1e16]")
        self.ve.br_in['column_height'] = float(a)

    @property
    def column_diameter(self):
        return self.ve.br_in['column_diameter']

    @column_diameter.setter
    def column_diameter(self, a):
        if not 1 <= a <= 6:
            raise ValueError(f"Value {a} is outside allowed interval [1, 1e16]")
        self.ve.br_in['column_diameter'] = float(a)

    @property
    def bubble_diameter(self):
        return self.ve.br_in['bubble_diameter']

    @bubble_diameter.setter
    def bubble_diameter(self, a):
        if not 0.003 <= a <= 0.008:
            raise ValueError(f"Value {a} is outside allowed interval [1, 1e16]")
        self.ve.br_in['bubble_diameter'] = float(a)

    @property
    def t_final(self):
        return self.ve.br_in['t_final']

    @t_final.setter
    def t_final(self, a):
        if not 1 <= a <= 1e16:
            raise ValueError(f"Value {a} is outside allowed interval [1, 1e16]")
        self.ve.br_in['t_final'] = float(a)

    @property
    def model_type(self):
        return self.ve.br_in['model_type']

    @model_type.setter
    def model_type(self, a):
        if not a in ['CFD Simulation', "CFD Surrogate"]:
            raise ValueError("Invalid value. Allowed options: 'CFD Simulation', 'CFD Surrogate'")
        self.ve.br_in['model_type'] = a
        self.select_run_function()
    ##############################################
    #
    ##############################################

    def select_run_function(self):
        # selected enzymatic hydrolysis model
        if self.model_type == 'CFD Simulation':
            assert self.hpc_run, f'Cannot run bioreactor without HPC resources. \n {os.getcwd()}'
            self.run = self.run_biorector_cfd_simulation
        elif self.model_type == "CFD Surrogate":
            self.run = self.run_biorector_cfd_surrogate
            br_surrogate_path = os.path.join(self.br_module_path,'surrogate_model')
            if not br_surrogate_path in sys.path:
                sys.path.append(br_surrogate_path)

    def run_biorector_cfd_simulation(self, verbose=True):

        print('\nRunning Bioreactor')
        os.chdir('bioreactor/bubble_column/')
        # Make changes to the fvOptions file based on replacement options
        fvOptions = {}
        fvOptions['rho_g'] = self.ve.eh_in['rho_g']
        fvOptions['rho_x'] = self.ve.eh_in['rho_x']
        fvOptions['rho_f'] = self.ve.eh_in['rho_f']
        write_file_with_replacements('constant/fvOptions', fvOptions)
        
        # Make changes to the controlDict file based on replacement options
        controlDict = {}
        controlDict['endTime'] = self.t_final
        write_file_with_replacements('system/controlDict', controlDict)

        # Run the bioreactor model
        # call function to update ovOptions # fvOptions?
        command = "sbatch ofoamjob"
        subprocess.run(command.split())
        self.ve.br_out = {'OUR': np.nan}
        os.chdir(self.notebookDir)
        print('Finished Bioreactor')

        if check_dict_for_nans(self.ve.br_out):
            return True
        return False


    def run_biorector_cfd_surrogate(self, verbose=True):
        
        print('\nRunning Bioreactor')
        
        if np.isnan(self.ve.eh_out['rho_g']):
            # TODO: Why we don't check it in CFD run_biorector_cfd_simulation? OD
            print('Waiting for EH CFD results.')
        else:
            from bcolumn_surrogate import run_br_surrogate
            self.ve.br_out = run_br_surrogate(self.ve)
            
            print('Finished Bioreactor')
            if check_dict_for_nans(self.ve.br_out):
                return True
            return False

def run_script(filename, *args, verbose=True):
    """ Execute the contents of a file.

    This function will attempt to execute the contents of a file specified
    with ``filename`` using the Python ``exec`` function.  No error checking
    is performed on the source file to be executed.

    Args:
        filename (str):
            The filename to execute line by line.
        *args:
            Variable length argument list to be made
            available to the executed file via ``sys.argv[..]``.
        verbose (bool, optional):
            Flag to display the printed outputs
            from the executed file, defaults to ``True``.

    Returns:
        None

    """

    sys.argv = [filename]
    sys.argv.extend(args)
    exec_file = open(filename, 'r')

    if verbose:
        # Execute the file as usual
        exec(exec_file.read(), globals())

    else:
        # Execute the file, redirecting all output to devnull
        # This suppresses any print statements within `filename`
        with open(os.devnull, 'w') as fp:
            with contextlib.redirect_stdout(fp):
                exec(exec_file.read(), globals())

    exec_file.close()