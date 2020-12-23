import subprocess
import yaml
import sys

def get_host_computer():

    hpc_run = False

    filename = 'hostname_info.txt'
    bash_command = 'hostnamectl'

    output_file = open(filename, 'w')
    try:
        process = subprocess.run(bash_command.split(), stdout=output_file)
    except:
        output_file.write('No hostname command found.')
    output_file.close()


    input_file = open(filename, 'r')
    for line in input_file:
        if 'computer-server' in line:
            hpc_run = True
            break
    input_file.close()

    if not hpc_run:
        print('It looks like you\'re running this notebook on a laptop.')
        print('Some operations requiring HPC resources will be disabled.')

    return hpc_run


def print_dict(dict_to_print, indent=0):

    for key, value in dict_to_print.items():
        if type(value) is dict:
            for j in range(indent):
                print('  ', end='')
            print('["%s"]' % (key))
            indent += 1
            print_dict(value, indent)
            indent -= 1

        else:
            for j in range(indent):
                print('  ', end='')
            print('["%s"] = ' % (key), value)

def dict_to_yaml(dictionary_to_write, yaml_filename, merge_with_existing=False, verbose=False):
    if type(dictionary_to_write) is not dict and type(dictionary_to_write) is not list:
        error_string = ('Input is not a dictionary. \nArguments should be (dict_to_write, yaml_filename), \nare they switched?')
        raise ValueError(error_string)

    if type(dictionary_to_write) is list:
        dict_temp = {}

        for d in dictionary_to_write:
            dict_temp.update(d)

        dictionary_to_write = dict_temp

    if merge_with_existing:
        existing_dict = yaml_to_dict(yaml_filename)
        existing_dict.update(dictionary_to_write)
        dictionary_to_write = existing_dict

    if verbose:
        print('Writing the Dictionary:')
        print_dict(dictionary_to_write)
        print('To the File: %s\n' % (yaml_filename))

    with open(yaml_filename, 'w') as fp:
        yaml.dump(dictionary_to_write, fp, sort_keys=False)


def yaml_to_dict(yaml_filename, verbose=False):
    with open(yaml_filename) as fp:
        output_dictionary = yaml.load(fp, Loader = yaml.FullLoader)

    if verbose:
        print('Read the Dictionary:')
        print_dict(output_dictionary)
        print('From the File: %s\n' % (yaml_filename))

    return output_dictionary


def run_script(filename, *args):
    """ Execute the contents of a file (`filename`). Optional arguments may be provided. """
    sys.argv = [filename]
    sys.argv.extend(args)
    file = open(filename)
    exec(file.read(), globals())
    file.close()
    return
