import subprocess
import yaml

def get_host_computer():

    found_hostname_command = True

    try:
        process = subprocess.Popen(['hostnamectl'], stdout=subprocess.PIPE)
    except:
        process = subprocess.Popen(['echo'], stdout=subprocess.PIPE)
        found_hostname_command = False


    hpc_run = False

    if found_hostname_command:
        command_line_output, error = process.communicate()

        hpc_test = [line for line in command_line_output if 'computer-server' in line]

        if len(hpc_test) > 0:
            hpc_run = True


    if not hpc_run:
        print('It looks like you\'re running this notebook on a laptop.')
        print('Some operations requiring HPC resources will be disabled.')

    return hpc_run


def dict_to_yaml(dictionary_to_write, yaml_filename, verbose=True):
    if type(dictionary_to_write) is not dict:
        error_string = ('Input is not a dictionary. \nArguments should be (dict_to_write, yaml_filename), \nare they switched?')
        raise ValueError(error_string)

    if verbose:
        print('Writing the Dictionary:')
        for k, v in dictionary_to_write.items():
            print('    ["%s"] = %f' % (k, v))
        print('To the File: %s\n' % (yaml_filename))

    with open(yaml_filename, 'w') as fp:
        yaml.dump(dictionary_to_write, fp)


def yaml_to_dict(yaml_filename, verbose=True):
    with open(yaml_filename) as fp:
        output_dictionary = yaml.load(fp, Loader = yaml.FullLoader)

    if verbose:
        print('Read the Dictionary:')
        for k, v in output_dictionary.items():
            print('    ["%s"] = %f' % (k, v))
        print('From the File: %s\n' % (yaml_filename))

    return output_dictionary
