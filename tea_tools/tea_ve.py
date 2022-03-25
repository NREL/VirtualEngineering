'''
This script calls the an Aspen sugar model and excel calculation
as part of the Virtual Engineering workflow.
'''

import sys
import os
import numpy as np

from classes import Aspen, Excel
from vebio.Utilities import dict_to_yaml, yaml_to_dict


def _init_tea_replacements(enzyme_loading=10.0, glucose_conv=0.9630, xylose_conv=0.9881, arabinose_conv=0.9881):
    tea_replacements = {}

    max_digits = 4

    tea_replacements['enzyme_loading'] = [os.path.join('Data', 'Flowsheeting Options', 'Calculator', 'A400-ENZ',
                                          'Input', 'FORTRAN_EXEC', '#22'),
                                          np.round(enzyme_loading, max_digits), True]

    tea_replacements['glucose_conv'] = [os.path.join('Data', 'Blocks', 'A300', 'Data', 'Blocks', 'CEH', 'Data',
                                        'Flowsheeting Options', 'Calculator', 'CEH', 'Input', 'FORTRAN_EXEC', '#11'),
                                        np.round(glucose_conv, max_digits), True]

    tea_replacements['xylose_conv'] = [os.path.join('Data', 'Blocks', 'A300', 'Data', 'Blocks', 'CEH', 'Data',
                                       'Flowsheeting Options', 'Calculator', 'CEH', 'Input', 'FORTRAN_EXEC', '#12'),
                                       np.round(xylose_conv, max_digits), True]

    tea_replacements['arabinose_conv'] = [os.path.join('Data', 'Blocks', 'A300', 'Data', 'Blocks', 'CEH', 'Data',
                                          'Flowsheeting Options', 'Calculator', 'CEH', 'Input', 'FORTRAN_EXEC', '#13'),
                                          np.round(arabinose_conv, max_digits), True]

    return tea_replacements

def run_aspen_model_ve(aspenFile, ve_params, outDir):

    try:
        # ================================================================
        # Create Aspen Plus communicator
        # ================================================================
        print('Opening Aspen Plus model... ')
        aspenModel = Aspen(aspenFile)
        print('Success!')

        # ================================================================
        # Write backup file with value updates/replacements
        # ================================================================
        # aspenPath: str, path in ASPEN tree
        # value: float or str, value to set
        # ifFortran: bool, whether it is a Fortran variable
        glucose_conv = ve_params['CEH_output']['System level'] ['Glucan conversion yield']
        xylose_conv = ve_params['CEH_output']['System level'] ['Xylan conversion yield']
        tea_replacements = _init_tea_replacements(glucose_conv=glucose_conv, xylose_conv=xylose_conv)

        print('Changing values in model backup definition tree... ')
        for key, val in tea_replacements.items():
            path_to_value = val[0]
            value = val[1]
            ifFortran = val[2]
            aspenModel.set_value(path_to_value, value, ifFortran, short_name=key)
        print('Success!')

        # ================================================================
        # Run the Aspen Plus model
        # ================================================================
        print('Running Aspen Plus model... ')
        aspenModel.run_model()
        print('Success!')

        # ================================================================
        # Save current model state
        # ================================================================
        print('Saving current model definition... ')
        tmpFile = os.path.join(outDir, 'CEH_ve.bkp')
        aspenModel.save_model(tmpFile)
        print('Success!')

    finally:
        aspenModel.close()

    return tmpFile


def run_excel_calc_ve(excelFile, tmpFile, ve_params):

    # ================================================================
    # Create Excel communicator and run calculator
    # ================================================================
    try:
        tmpFile = os.path.abspath(tmpFile)
    
        print('Opening Excel calculator... ')
        excelCalculator = Excel(excelFile)
        excelCalculator.load_aspenModel(tmpFile)
        print('Success!')

        power = excelCalculator.get_cell('OPEX', 'E40')
        print(f'Old Power: {power}')
        delta_power = ve_params['CEH_output']['System level']['Total membrane loop power consumption (kW)']
        print(f'Delta Power: {delta_power}')
        power += delta_power

        excelCalculator.set_cell(power, 'OPEX', 'E40')
        power = excelCalculator.get_cell('OPEX', 'E40')
        print(f'New Power: {power}')

        print('Running Excel analysis... ')
        excelCalculator.run_macro('solvedcfror')
        print('Success!')
        
        mssp = excelCalculator.get_cell('DCFROR', 'B36')

    finally:
        excelCalculator.close()

    return mssp
    

def run_tea_ve(*args):

    if len(sys.argv) > 1:
        params_filename = sys.argv[1]
        ve_params = yaml_to_dict(params_filename)
    else:
        raise Exception("VE parameters filename not provided")

    # Set the input and output files/directories
    aspenFile = os.path.abspath(ve_params['tea_input']['aspen_filename'])
    excelFile = os.path.abspath(ve_params['tea_input']['excel_filename'])

    outDir = 'output'
    os.makedirs(outDir, exist_ok = True)

    tmpFile = run_aspen_model_ve(aspenFile, ve_params, outDir)

    mssp = run_excel_calc_ve(excelFile, tmpFile, ve_params)

    print(f'Selling price: {mssp:.6f}')


run_tea_ve(sys.argv)
