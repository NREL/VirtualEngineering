import os
import yaml

def get_host_computer():
    """ Check if environment is running on the HPC.

    This function queries an environment variable, ``NREL_CLUSTER`` and if the pattern
    "eagle" is found, sets the variable ``hpc_run`` to True to 
    allow the submission CFD jobs from the Notebook.  Note that the particular
    pattern "eagle" is found to correctly identify all instances
    of the Eagle HPC but may need to be updated to catch other
    HPC systems (e.g., Vermillion or Swift).
    
    Args:
        None

    Returns:
        bool:
            True if HPC resources identified, False otherwise. 

    """

    nrel_cluster = os.environ.get('NREL_CLUSTER')

    if nrel_cluster == 'eagle':
        hpc_run = True

    else:
        hpc_run = False

        print('It looks like you\'re running this notebook on a laptop.')
        print('Operations requiring HPC resources will be disabled.')

    return hpc_run


def print_dict(dict_to_print, indent=0):
    """ Neatly print a dictionary with indentation for nesting.

    This function prints a dictionary with items visually indented according to 
    their nesting.  Multiple nesting and sub-nesting levels are supported.
    
    Args:
        dict_to_print (dict):
            The dictionary to print. The dictionary entries
            can be sub-dictionaries which creates the nested structure.
        indent (int, optional):
            The indentation level to use at the root level of
            the dictionary.  If not included, the default indentation of 0 is used. 
            Non-zero values of indent effectively shift the output to the right.

    Returns:
        None

    """

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
    """ Write an input dictionary to a YAML file.

    This function takes an input dictionary and writes it directly to
    a YAML file, preserving any nesting and keyname relationships.  For
    example an input dictionary::

        input_dict = {'parent_1': {'a': 2, 'b': 3},
                      'parent_2': {'a': 20, 'b': 30}}

    will be translated into a YAML file with the same hierarchy::

        parent_1:
          a: 2
          b: 3
        parent_2:
          a: 20
          b: 30

    This function also allows for the converted dictionary to be appended
    to an existing YAML file.
    
    Args:
        dictionary_to_write (dict, list(dict)):
            The dictionary or a list of dictionaries to save as a YAML file.
        yaml_filename (str):
            The filename to use when saving the YAML file.
        merge_with_existing (bool, optional):
            A flag to specify that the ``yaml_filename`` already exists 
            and the new values help in ``dictionary_to_write`` should be
            appended to this existing file rather than overwriting it.
            Defaults to ``False``.
        verbose (bool, optional):
            A flag to print additional output, including the contents of the 
            dictionary, after the conversion process. Defaults to ``False``

    Returns:
        None

    """

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
    """ Read a YAML file and store as a dictionary.

    This function reads an input YAML file and stores the contents
    as a Python dictionary, where all necessary relationships and 
    keynames are preserved as shown in the YAML file.
    
    Args:
        yaml_filename (str):
            The YAML filename to read the data from.
        verbose (bool, optional):
            A flag to print additional output, including the contents of the 
            converted dictionary, after the conversion process. Defaults to ``False``

    Returns:
        None

    """

    with open(yaml_filename) as fp:
        output_dictionary = yaml.load(fp, Loader = yaml.FullLoader)

    if verbose:
        print('Read the Dictionary:')
        print_dict(output_dictionary)
        print('From the File: %s\n' % (yaml_filename))

    return output_dictionary
