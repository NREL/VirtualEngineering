import sys
import os
import shutil
import contextlib
import subprocess
import glob

from scipy.interpolate import interp1d
import numpy as np

from vebio.FileModifiers import write_file_with_replacements
from vebio.Utilities import yaml_to_dict, dict_to_yaml

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

    method = 1

    if method == 1:
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

    elif method == 2:
        # Begin the command with `python <filename>`
        command = f'python {filename}'

        # Extend with any command-line arguments
        for arg in args:
            command += f' {arg}'

        # Execute this command using subprocess
        command_output = subprocess.run(command.split(), capture_output=True, text=True)

        if verbose:
            print(command_output.stdout)

    return

def run_pretreatment(notebookDir, params_filename, fs_options, pt_options, verbose=True):
    """ Run the pretreatment operation.

    This function runs the pretreatment unit model specified in ``ptrun.py``.
    Since this is the first unit operation in the overall conversion process,
    the feedstock properties are integrated during this step.

    Through the ``fs_options`` widgets, the user controls the following
    values:

        * The initial fraction of solids due to xylan (X_X)
        * The initial fraction of solids due to glucan (X_G)
        * The initial porous fraction of the biomass particles

    Through the ``pt_options`` widgets, the user controls the following
    values:

        * Acid Loading (float)
        * Steam Temperature (float)
        * Initial FIS_0 (float)
        * Final Time (float)
        * Show plots (bool)

    Args:
        notebookDir (str):
            The path to the Jupyter Notebook, used to specify the location
            of the input file and reset the working directory after this operation
            is finished.

        params_filename (str):
            The filename for the parameters yaml file including
            extension, e.g., ``'virteng_params.yaml'``

        fs_options (WidgetCollection):
            A ``WidgetCollection`` object containing all of widgets used
            to solicit user input for feedstock properties.

        pt_options (WidgetCollection):
            A ``WidgetCollection`` object containing all of widgets used
            to solicit user input for pretreatment properties.

        verbose (bool, optional):
            Option to show print messages from executed file, default True.

    Returns:
        None

    """

    print('Running Pretreatment Model')
    
    # Export the feedstock and pretreatment options to a global yaml file
    fs_dict = fs_options.export_widgets_to_dict(parent_name='feedstock')
    pt_dict = pt_options.export_widgets_to_dict(parent_name='pretreatment_input')

    # Obtain steam concentration from lookup table and add to dictionary
    steam_data = np.genfromtxt('pretreatment_model/lookup_tables/sat_steam_table.csv', delimiter=',', skip_header=1)

    # build interpolator interp_steam = interp.interp1d(temp_in_K, dens_in_kg/m3)
    interp_steam = interp1d(steam_data[:, 2], steam_data[:, 4])
    dens = interp_steam(pt_dict['pretreatment_input']['steam_temperature'])

    # Convert to mol/ml => density in g/L / molecular weight / 1000.0
    mol_per_ml = float(dens/18.01528/1000.0)
    pt_dict['pretreatment_input']['bulk_steam_conc'] = mol_per_ml

    dict_to_yaml([fs_dict, pt_dict], params_filename)

    # Move into the pretreatment directory
    test_folder_path = os.path.join(notebookDir, 'pretreatment_model/test/')
    sys.path.append(test_folder_path)
    os.chdir(test_folder_path)
    
    # See if the pretreatment module exists, if not, we need to build it
    try:
        import pt
    except:
        print('Could not load PT module, building module from source.')
        print('(This will only happen the first time the notebook is run.)')
        build_folder_path = os.path.join(notebookDir, 'pretreatment_model/bld/')
        command = f'make -C {build_folder_path} ptpython'
        subprocess.run(command.split())
        print('Finished building PT module.')

        print('Copying modules to test directory.')
        files_to_copy = glob.glob(f'{build_folder_path}*.so')
        print(f'Found {len(files_to_copy)} files to copy:')
        for f in files_to_copy:
            print(f'| Copying {f} to directory {test_folder_path}')
            shutil.copy(f, test_folder_path)
        print('Finished copying files.')

    os.chdir(test_folder_path)

    # clear out old data files (`postprocess.py` will pick up longer-run stale data files)
    outfiles = glob.glob("out*.dat")
    for outfile in outfiles:
        os.remove(outfile)
    # Run pretreatment code specifying location of input file
    path_to_input_file = os.path.join(notebookDir, params_filename)
    for k, p in enumerate(sys.path):
        print(k, p)
    run_script("ptrun.py", path_to_input_file, verbose=verbose)
    # unwinding the below because a fix to `f2pymain.f90` now allows rerunning
    # `ptrun.py`; not sure if capturing the output is still wanted, though; JJS
    # 1/13/21
    #pt_run_command = 'python ptrun.py %s' % (path_to_input_file)
    #pt_cli = subprocess.run(pt_run_command.split(), capture_output=True, text=True)
    #print(pt_cli.stdout[-1394:])

    if pt_options.show_plots.value:
        run_script("postprocess.py", "out_*.dat", "exptdata_150C_1acid.dat", verbose=verbose)

    os.chdir(notebookDir)
    print('\nFinished Pretreatment')

    
def run_enzymatic_hydrolysis(notebookDir, params_filename, eh_options, hpc_run,
                             verbose=True):
    """ Run the enzymatic hydrolysis operation.

    This function runs the enzymatic hydrolysis unit operation. Three 
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

    Args:
        notebookDir (str):
            The path to the Jupyter Notebook, used to specify the location
            of the input file and reset the working directory after this operation
            is finished.

        params_filename (str):
            The filename for the parameters yaml file including
            extension, e.g., ``'virteng_params.yaml'``

        eh_options (WidgetCollection):
            A ``WidgetCollection`` object containing all of widgets used
            to solicit user input for enzymatic hydrolysis properties.

        hpc_run (bool):
            A flag indicating whether or not the Notebook is being
            run on HPC resources, enable CFD only if True.

        verbose (bool, optional):
            Option to show print messages from executed file, default True.

    Returns:
        None

    """

    print('\nRunning Enzymatic Hydrolysis Model')

    # Export the enzymatic hydrolysis options to a global yaml file
    eh_dict = eh_options.export_widgets_to_dict(parent_name='enzymatic_input')
    dict_to_yaml(eh_dict, params_filename, merge_with_existing=True)
    
    # Run the selected enzymatic hydrolysis model
    if eh_options.model_type.value == 'CFD Simulation':
        
        # Export the current state to a dictionary
        ve_params = yaml_to_dict(params_filename)
        os.chdir('EH_OpenFOAM/tests/RushtonReact/')
                
        # Prepare input values for EH CFD operation
        globalVars = {}

        globalVars['fis0'] = ve_params['enzymatic_input']['fis_0']
        globalVars['xG0'] = ve_params['pretreatment_output']['X_G']
        globalVars['xX0'] = ve_params['pretreatment_output']['X_X']
        globalVars['XL0'] = 1.0 - globalVars['xG0'] - globalVars['xX0']
        globalVars['yF0'] = 0.2 + 0.6*ve_params['pretreatment_output']['conv']
        globalVars['lmbdE'] = ve_params['enzymatic_input']['lambda_e']
        globalVars['rhog0'] = 0.0
        dilution_factor = ve_params['enzymatic_input']['fis_0']/ve_params['pretreatment_output']['fis_0']
        globalVars['rhox0'] = ve_params['pretreatment_output']['rho_x']*dilution_factor
        globalVars['rhosl0'] = 0.0

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
        fintime = fluid_steadystate_time + (ve_params['enzymatic_input']['t_final']/reaction_update_time + 1.0)*fluid_update_time
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


        if hpc_run:
            # command = "srun hostname"
            # host_list = subprocess.run(command.split(), capture_output=True).stdout.decode()
            # num_nodes = len(host_list)
            # max_cores = int(36*num_nodes)

            username = os.environ['USER']
            jobname = 'eh_cfd'

            command = 'squeue -u %s -t R,PD -n %s' % (username, jobname)
            out = subprocess.run(command.split(), capture_output=True, text=True)

            output_dict = {'enzymatic_output': {}}
            output_dict['enzymatic_output']['rho_g'] = np.nan
            output_dict['enzymatic_output']['rho_x'] = np.nan
            output_dict['enzymatic_output']['rho_sl'] = np.nan
            output_dict['enzymatic_output']['rho_f'] = np.nan

            if username in out.stdout:
                # Job is running, do nothing
                print('EH CFD job is already queued.')
                print(out.stdout)
                job_id = out.stdout.strip().split('\n')[-1].split()[0]

                if eh_options.use_previous_output.value:
                    print('Using outputs from most recent finished simulation.')
                    integrated_quantities = np.genfromtxt('old_integrated_quantities.dat') # mol/L
                    output_dict = {'enzymatic_output': {}}
                    output_dict['enzymatic_output']['rho_g'] = float(integrated_quantities[-1, -3])
                    output_dict['enzymatic_output']['rho_x'] = float(integrated_quantities[-1, -2])
                    output_dict['enzymatic_output']['rho_sl'] = float(integrated_quantities[-1, -1])
                    output_dict['enzymatic_output']['rho_f'] = float(ve_params['pretreatment_output']['rho_f']*dilution_factor)
                    print('Success.')

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

        else:
            print('Cannot run EH_CFD without HPC resources.')
            print(os.getcwd())
            
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

        
        # output_dict = {'enzymatic_output': {}}
        # output_dict['enzymatic_output']['rho_g'] = integrated_quantities[-1, -3]
        # output_dict['enzymatic_output']['rho_x'] = integrated_quantities[-1, -2]
        # output_dict['enzymatic_output']['rho_sl'] = integrated_quantities[-1, -1]
        # output_dict['enzymatic_output']['rho_f'] = ve_params['pretreatment_output']['rho_f']*dilution_factor
        
        os.chdir(notebookDir)

        dict_to_yaml([ve_params, output_dict], params_filename)
        
    else:
        
        path_to_input_file = os.path.join(notebookDir, params_filename)

        if eh_options.model_type.value == 'CFD Surrogate':
            os.chdir('EH_OpenFOAM/EH_surrogate/')
            run_script("EH_surrogate.py", path_to_input_file, verbose=verbose)

        else:
            os.chdir('two_phase_batch_model/')
            # Commenting out cellulose-only two-phase model to use lignocellulose
            # model, just in case we want to switch back or make both an
            # option. The lignocellulose model is superior.
            #run_script("two_phase_batch_model.py", path_to_input_file, verbose=verbose)
            run_script("driver_batch_lignocell_EH_VE.py", path_to_input_file, verbose=verbose)
        
        os.chdir(notebookDir)

    print('\nFinished Enzymatic Hydrolysis')

    
def run_bioreactor(notebookDir, params_filename, br_options, hpc_run, verbose=True):
    """ Run the aerobic bioreaction operation.

    This function runs the aerobic bioreaction operation using the
    user-specified properties.  Two distinct models exist: (1) a 
    pre-trained surrogate model informed from CFD runs and (2) the
    full CFD simulation itself where option (2) is accessible only
    with ``hpc_run=True``.  The default option is the surrogate model.

    Through the ``br_options`` widgets, the user controls the following
    values:

        * Model Type
        * Final Time (float)

    Args:
        notebookDir (str):
            The path to the Jupyter Notebook, used to specify the location
            of the input file and reset the working directory after this operation
            is finished.

        params_filename (str):
            The filename for the parameters yaml file including
            extension, e.g., ``'virteng_params.yaml'``

        br_options (WidgetCollection):
            A ``WidgetCollection`` object containing all of widgets used
            to solicit user input for bioreaction properties.

        hpc_run (bool):
            A flag indicating whether or not the Notebook is being
            run on HPC resources, enable CFD only if True.

        verbose (bool, optional):
            Option to show print messages from executed file, default True.

    Returns:
        None

    """

    print('\nRunning Bioreactor Model')

    # Export the bioreactor options to a global yaml file
    br_dict = br_options.export_widgets_to_dict(parent_name='bioreactor_input')
    dict_to_yaml(br_dict, params_filename, merge_with_existing=True)

    # Convert the current parameters file to a dictionary
    ve_params = yaml_to_dict(params_filename)

    # Run the selected CFD or surrogate model
    if br_options.model_type.value == 'CFD Simulation':

        os.chdir('bioreactor/bubble_column/')

        # Make changes to the fvOptions file based on replacement options
        fvOptions = {}

        fvOptions['rho_g'] = ve_params['enzymatic_output']['rho_g']
        fvOptions['rho_x'] = ve_params['enzymatic_output']['rho_x']
        fvOptions['rho_f'] = ve_params['enzymatic_output']['rho_f']

        write_file_with_replacements('constant/fvOptions', fvOptions)
        
        # Make changes to the controlDict file based on replacement options
        controlDict = {}
        controlDict['endTime'] = ve_params['bioreactor_input']['t_final']
        
        write_file_with_replacements('system/controlDict', controlDict)

        # Run the bioreactor model
        if hpc_run:
            # call function to update ovOptions # fvOptions?
            command = "sbatch ofoamjob"
            subprocess.run(command.split())
        else:
            print('Cannot run bioreactor without HPC resources.')
            print(os.getcwd())
            
        output_dict = {'bioreactor_output': {}}
        output_dict['bioreactor_output']['placeholder'] = 123

        os.chdir(notebookDir)

        dict_to_yaml([ve_params, output_dict], params_filename)
        
    else:
        if np.isnan(ve_params['enzymatic_output']['rho_g']):
            print('Waiting for EH CFD results.')
        else:
            os.chdir('bioreactor/bubble_column/surrogate_model')
            path_to_input_file = os.path.join(notebookDir, params_filename)
            run_script("bcolumn_surrogate.py", path_to_input_file, verbose=verbose)
            os.chdir(notebookDir)

    print('\nFinished Bioreactor')
