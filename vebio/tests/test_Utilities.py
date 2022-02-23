import pytest
import os
from contextlib import redirect_stdout
from io import StringIO
import yaml

from vebio.Utilities import get_host_computer, print_dict, dict_to_yaml, yaml_to_dict

test_yaml_filename = 'temp.yaml'

def test_get_host_computer():
    # Test with the current environment variable
    hpc_run = get_host_computer()
    assert hpc_run == False

    # Set the environment variable by hand
    os.environ['NREL_CLUSTER'] = 'eagle'
    hpc_run = get_host_computer()
    assert hpc_run == True

    # Clean up the manually-set variable
    del os.environ['NREL_CLUSTER']

def test_print_dict():

    input_dict = {'A': 100, 'B': {'x': 10, 'y': 20, 'z': {'m': 1, 'n': 2}}, 'C': 200}

    fp = StringIO()

    with redirect_stdout(fp):
        print_dict(input_dict)

    output = fp.getvalue().split('\n')

    truth_value = ['["A"] =  100',
                   '["B"]',
                   '  ["x"] =  10',
                   '  ["y"] =  20',
                   '  ["z"]',
                   '    ["m"] =  1',
                   '    ["n"] =  2',
                   '["C"] =  200']
        
    for k, j in zip(output, truth_value):
        assert k == j

def test_dict_to_yaml():
    
    input_dict = {'model_type_1': 'CFD Surrogate', 'initial_porosity_1': 0.234, 'show_plots_1': False}

    dict_to_yaml(input_dict, test_yaml_filename)

    with open(test_yaml_filename) as fp:
        output_dict = yaml.load(fp, Loader = yaml.FullLoader)

    for key, val in input_dict.items():
        assert key in output_dict
        assert input_dict[key] == output_dict[key]

def test_dict_to_yaml_multi_merge():
    
    input_dict_a = {'model_type_2': 'Batch Model', 'initial_porosity_2': 0.345, 'show_plots_2': True}
    input_dict_b = {'model_type_3': 'CFD',         'initial_porosity_3': 0.456, 'show_plots_3': False}

    dict_to_yaml([input_dict_a, input_dict_b], test_yaml_filename, merge_with_existing=True)

    with open(test_yaml_filename) as fp:
        output_dict = yaml.load(fp, Loader = yaml.FullLoader)

    print(output_dict)

    assert len(output_dict) == 9

    for key, val in input_dict_a.items():
        assert key in output_dict
        assert val == output_dict[key]

    for key, val in input_dict_b.items():
        assert key in output_dict
        assert val == output_dict[key]

def test_yaml_to_dict():
    
    output_dict = yaml_to_dict(test_yaml_filename)

    truth_dict = {'model_type_1': 'CFD Surrogate', 'initial_porosity_1': 0.234, 'show_plots_1': False,
                  'model_type_2': 'Batch Model',   'initial_porosity_2': 0.345, 'show_plots_2': True,
                  'model_type_3': 'CFD',           'initial_porosity_3': 0.456, 'show_plots_3': False}

    for key, val in truth_dict.items():
        assert key in output_dict
        assert val == output_dict[key]

    os.remove(test_yaml_filename)
