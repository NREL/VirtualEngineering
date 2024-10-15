import pytest
import os

from virteng.FileModifiers import write_file_with_replacements

@pytest.fixture()
def build_test_file():

    filename = 'ci_temp.txt'

    with open(filename, 'w') as fp:
        fp.write('A = 1.0;\n')
        fp.write('B: 1.0 \n')
        fp.write('C    1.0\n')
        fp.write(' D 1.0\n')
        fp.write(' DD 1.0 \n')

    return filename

@pytest.fixture()
def build_truth_file():

    filename = 'ci_temp_truth.txt'

    with open(filename, 'w') as fp:
        fp.write('A = 2.2;\n')
        fp.write('B: 2.2\n')
        fp.write('C 2.2\n')
        fp.write('D 2.2\n')
        fp.write(' DD 1.0 \n')

    return filename

@pytest.fixture()
def build_test_dict():

    replacement_dict = {}

    replacement_dict['A'] = 2.2
    replacement_dict['B'] = 2.2 
    replacement_dict['C'] = 2.2
    replacement_dict['D'] = 2.2

    return replacement_dict

@pytest.mark.unit
def test_replacement(build_test_file, build_test_dict, build_truth_file):

    test_file = build_test_file
    truth_file = build_truth_file
    test_dict = build_test_dict

    write_file_with_replacements(test_file, test_dict)

    with open(test_file, 'r') as f_test, open(truth_file, 'r') as f_truth:
        test_data = f_test.readlines()
        truth_data = f_truth.readlines()

    assert len(test_data) == len(truth_data)

    for test_val, truth_val in zip(test_data, truth_data):
        assert test_val == truth_val

    # Clean up
    os.remove(test_file)
    os.remove('bkup_'+test_file)
    os.remove(truth_file)

