'''
This script calls the an Aspen sugar model and excel calculation
as part of the Virtual Engineering workflow.
'''

import sys
import os
import numpy as np

from classes import Aspen, Excel
from vebio.Utilities import dict_to_yaml, yaml_to_dict


def _init_tea_replacements(enzyme_loading=12.0, glucose_conv=0.9630, xylose_conv=0.9881, recirc_ratio=6.5):
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

    tea_replacements['recirc_ratio'] = [os.path.join('Data', 'Blocks', 'A300', 'Data', 'Blocks', 'CEH', 'Data',
                                       'Blocks', 'B19', 'Input', 'FACTOR'),
                                       np.round(recirc_ratio, max_digits), False]

    # tea_replacements['delta_p'] = [os.path.join('Data', 'Blocks', 'A300', 'Data', 'Blocks', 'CEH', 'Data',
    #                                    'Blocks', 'B18', 'Input', 'DELP'),
    #                                    np.round(delta_p, max_digits), False]

    # tea_replacements['arabinose_conv'] = [os.path.join('Data', 'Blocks', 'A300', 'Data', 'Blocks', 'CEH', 'Data',
    #                                       'Flowsheeting Options', 'Calculator', 'CEH', 'Input', 'FORTRAN_EXEC', '#13'),
    #                                       np.round(arabinose_conv, max_digits), True]

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
        recirc_ratio = ve_params['CEH_output']['System level']['Overall retentate/permeate ratio']
        tea_replacements = _init_tea_replacements(glucose_conv=glucose_conv, xylose_conv=xylose_conv, recirc_ratio=recirc_ratio)

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
        tmpFile = os.path.abspath(tmpFile)
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
    
        print('Opening Excel calculator... ')
        excelCalculator = Excel(excelFile)
        excelCalculator.load_aspenModel(tmpFile)
        print('Success!')

        size_1 = excelCalculator.get_cell('CAPEX', 'G70')
        size_2 = excelCalculator.get_cell('CAPEX', 'G76')
        size_3 = excelCalculator.get_cell('CAPEX', 'G82')
        mf_area = excelCalculator.get_cell('CAPEX', 'O91')

        size_1_new = ve_params['CEH_output']['CEH Reactor 1']['Reactor Size (kg)']
        size_2_new = ve_params['CEH_output']['CEH Reactor 2']['Reactor Size (kg)']
        size_3_new = ve_params['CEH_output']['CEH Reactor 3']['Reactor Size (kg)']
        mf_area_new = ve_params['CEH_output']['System level']['Total membrane area (m2)']

        print(f'size_1 = {size_1_new} kg (originally {size_1})')
        print(f'size_2 = {size_2_new} kg (originally {size_2})')
        print(f'size_3 = {size_3_new} kg (originally {size_3})')
        print(f'mf_area = {mf_area_new} m^2 (originally {mf_area})')

        excelCalculator.set_cell(size_1_new, 'CAPEX', 'G70')
        excelCalculator.set_cell(size_2_new, 'CAPEX', 'G76')
        excelCalculator.set_cell(size_3_new, 'CAPEX', 'G82')
        excelCalculator.set_cell(mf_area_new, 'CAPEX', 'O91')

        print('Running Excel analysis... ')
        excelCalculator.run_macro('solvedcfror')
        print('Success!')
        
        mssp = excelCalculator.get_cell('DCFROR', 'B36')

        output_dict = {}
        output_dict['tea_output'] = {}
        output_dict['tea_output']['mssp'] = float(mssp)

    finally:
        excelCalculator.close()

    return output_dict
    

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

    output_dict = run_excel_calc_ve(excelFile, tmpFile, ve_params)

    print(f'Selling price: {output_dict["tea_output"]["mssp"]:.6f}')

    dict_to_yaml([ve_params, output_dict], params_filename)

run_tea_ve(sys.argv)
