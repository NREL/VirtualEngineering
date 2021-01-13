import os
import subprocess

from vebio.FileModifiers import write_file_with_replacements
from vebio.Utilities import yaml_to_dict, dict_to_yaml, run_script

def run_pretreatment(notebookDir, params_filename, fs_options, pt_options):
    print('Running Pretreatment Model')
    
    # Export the feedstock and pretreatment options to a global yaml file
    fs_dict = fs_options.export_widgets_to_dict('feedstock')
    pt_dict = pt_options.export_widgets_to_dict('pretreatment_input')
    dict_to_yaml([fs_dict, pt_dict], params_filename)

    # Move into the pretreatment directory
    os.chdir('pretreatment_model/test/')
    
    # See if the pretreatment module exists, if not, we need to build it
    try:
        import pt
    except:
        print('Could not load PT module, building module from source.')
        print('(This will only happen the first time the notebook is run.)')
        os.chdir('../bld/')
        command = "sh build_first_time.sh"
        subprocess.run(command.split())
        os.chdir('../test/')
        print('Finished building PT module.')
              
    # Run pretreatment code specifying location of input file
    path_to_input_file = os.path.join(notebookDir, params_filename)
    # run_script("ptrun.py", path_to_input_file)
    pt_run_command = 'python ptrun.py %s' % (path_to_input_file)
    pt_cli = subprocess.run(pt_run_command.split(), capture_output=True, text=True)
    print(pt_cli.stdout[-1394:])

    if pt_options.show_plots.value:
        run_script("postprocess.py", "out_*.dat", "exptdata_150C_1acid.dat")

    os.chdir(notebookDir)
    print('\nFinished Pretreatment')

    
def run_enzymatic_hydrolysis(notebookDir, params_filename, eh_options, hpc_run):
    print('\nRunning Enzymatic Hydrolysis Model')

    # Export the enzymatic hydrolysis options to a global yaml file
    eh_dict = eh_options.export_widgets_to_dict('enzymatic_input')
    dict_to_yaml(eh_dict, params_filename, merge_with_existing=True)
    
    # Run the selected enzymatic hydrolysis model
    if eh_options.use_cfd.value:
        
        # Export the current state to a dictionary
        ve_params = yaml_to_dict(params_filename)
                
        # Prepare input values for EH CFD operation
        dilution_factor = ve_params['enzymatic_input']['fis_0']/ve_params['pretreatment_output']['fis_0']
        rho_x0 = ve_params['pretreatment_output']['rho_x']*dilution_factor
        rho_f0 = ve_params['pretreatment_output']['rho_f']*dilution_factor
 
        enzdata_replacements = {}
        enzdata_replacements['lmbde'] = ve_params['enzymatic_input']['lambda_e']
        enzdata_replacements['fis0'] = ve_params['enzymatic_input']['fis_0']
        enzdata_replacements['yF0'] = 0.2 + 0.6*ve_params['pretreatment_output']['conv']
        enzdata_replacements['fGl0'] = ve_params['pretreatment_output']['X_G']
        
        os.chdir('EH_CFD/')
        write_file_with_replacements('enzdata', enzdata_replacements)
        
        # Get dt_ss, dt_react, and dt_fr from enzdata to calculate FINTIME
        fp = open('enzdata', 'r')
        dt_ss = 0.0
        dt_react = 0.0
        dt_fr = 0.0
        for line in fp:
            if '#' not in line:
                if 'dt_ss' in line:
                    dt_ss = float(line.split('=')[1])
                elif 'dt_react' in line:
                    dt_react = float(line.split('=')[1])
                elif 'dt_fr' in line:
                    dt_fr = float(line.split('=')[1])
        fp.close()
        
        fintime = dt_ss + (ve_params['enzymatic_input']['t_final']/dt_react)*dt_fr

        paddle_rea_replacements = {}
        paddle_rea_replacements['FINTIME'] = '   %.5f     p010 FINTIME 480\n' % (fintime)
        write_file_with_replacements('paddle.rea', paddle_rea_replacements, full_overwrite=True)

        if hpc_run:
            command = "srun hostname"
            host_list = subprocess.run(command.split(), capture_output=True).stdout.decode()
            num_nodes = len(host_list)
            max_cores = int(36*num_nodes)

            command = "sh run.sh $s" % max_cores
            subprocess.run(command.split())
        else:
            print('Cannot run EH_CFD without HPC resources.')
            print('$ ./run.sh $max_cores') # not sure of the purpose of this line, JJS
            print(os.getcwd())
            
        # Prepare output values from EH CFD operations
        # FIXME: rho_g should be value taken from CFD output - should be good now, please check, JJS 9/15/20
        
        # per Hari's email, glucose concentration in fort.44 is mol/L
        c_g_output = np.genfromtxt('fort.44') # mol/L
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
        
        output_dict = {'enzymatic_output': {}}
        output_dict['enzymatic_output']['rho_g'] = rho_g_final
        output_dict['enzymatic_output']['rho_x'] = rho_x_final
        output_dict['enzymatic_output']['rho_f'] = rho_f_final
        
        os.chdir(notebookDir)

        dict_to_yaml([ve_params, output_dict], params_filename)
        
    else:
        os.chdir('two_phase_batch_model/')
        path_to_input_file = os.path.join(notebookDir, params_filename)
        # this works in notebooks but wouldn't elsehwere
        #%run two_phase_batch_model.py $path_to_input_file 
        run_script("two_phase_batch_model.py", path_to_input_file)
        os.chdir(notebookDir)

    print('\nFinished Enzymatic Hydrolysis')

    
def run_bioreactor(notebookDir, params_filename, br_options, hpc_run):
    print('\nRunning Bioreactor Model')

    # Export the bioreactor options to a global yaml file
    br_dict = br_options.export_widgets_to_dict('bioreactor_input')
    dict_to_yaml(br_dict, params_filename, merge_with_existing=True)

    # Run the selected enzymatic hydrolysis model
    if br_options.use_cfd.value:
        # Convert the current parameters file to a dictionary
        ve_params = yaml_to_dict(params_filename)        

        # Make changes to the fvOptions file based on replacement options
        fvOptions_replacements = {}
        for key, value in ve_params['enzymatic_output'].items():
            fvOptions_replacements['double %s' % (key)] = value

        os.chdir('bioreactor/bubble_column/constant/')
        write_file_with_replacements('fvOptions', fvOptions_replacements)
        os.chdir(notebookDir)
        
        # Make changes to the controlDict file based on replacement options
        controlDict_replacements = {}
        controlDict_replacements['endTime '] = ve_params['bioreactor_input']['t_final']
        
        os.chdir('bioreactor/bubble_column/system/')
        write_file_with_replacements('controlDict', controlDict_replacements)
        os.chdir(notebookDir)

        # Run the bioreactor model
        os.chdir('bioreactor/bubble_column/')
        if hpc_run:
            # call function to update ovOptions # fvOptions?
            command = "sbatch ofoamjob"
            subprocess.run(command.split())
        else:
            print('Cannot run bioreactor without HPC resources.')
            print('$ sbatch ofoamjob') # not sure of the purpose of this line, JJS
            print(os.getcwd())
            
        output_dict = {'bioreactor_output': {}}
        output_dict['bioreactor_output']['placeholder'] = 123

        os.chdir(notebookDir)

        dict_to_yaml([ve_params, output_dict], params_filename)
        
    else:
        os.chdir('bioreactor/bubble_column/surrogate_model')
        path_to_input_file = '%s/%s' % (notebookDir, params_filename)
        # %run bcolumn_surrogate.py $path_to_input_file
        run_script("bcolumn_surrogate.py", path_to_input_file)
        os.chdir(notebookDir)

    print('\nFinished Bioreactor')
