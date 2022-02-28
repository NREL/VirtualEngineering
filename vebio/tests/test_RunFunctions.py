import pytest
import os
from ipywidgets import *
from pathlib import Path

from vebio.RunFunctions import run_pretreatment, run_enzymatic_hydrolysis, run_bioreactor
from vebio.WidgetFunctions import WidgetCollection
from vebio.Utilities import yaml_to_dict


notebook_dir = os.getcwd()
params_filename = 'test_params.yaml'


@pytest.fixture()
def build_fs_options():
    fs_options = WidgetCollection()

    fs_options.xylan_solid_fraction = widgets.BoundedFloatText(value = 0.263)
    fs_options.glucan_solid_fraction = widgets.BoundedFloatText(value = 0.40)
    fs_options.initial_porosity = widgets.BoundedFloatText(value = 0.8)

    return fs_options

@pytest.fixture()
def build_pt_options():
    pt_options = WidgetCollection()

    pt_options.initial_acid_conc = widgets.BoundedFloatText(value = 0.0001)
    pt_options.steam_temperature = widgets.BoundedFloatText(value = 150.0)
    pt_options.steam_temperature.scaling_fn = lambda C : C + 273.15
    pt_options.initial_solid_fraction = widgets.BoundedFloatText(value = 0.745)
    pt_options.final_time = widgets.BoundedFloatText(value = 8.3)
    pt_options.final_time.scaling_fn = lambda s : 60.0 * s
    pt_options.show_plots = widgets.Checkbox(value = False)

    return pt_options

@pytest.fixture()
def build_eh_options():
    eh_options = WidgetCollection()

    eh_options.model_type = widgets.RadioButtons(
        options = ['Lignocellulose Model', 'CFD Surrogate', 'CFD Simulation'],value = 'CFD Surrogate')
    eh_options.lambda_e = widgets.BoundedFloatText(value = 30.0,)
    eh_options.lambda_e.scaling_fn = lambda e : 0.001 * e
    eh_options.fis_0 = widgets.BoundedFloatText(value = 0.05)
    eh_options.t_final = widgets.BoundedFloatText(value = 24.0)
    eh_options.show_plots = widgets.Checkbox(value = False)

    return eh_options

@pytest.fixture()
def build_br_options():
    br_options = WidgetCollection()

    br_options.model_type = widgets.RadioButtons(
        options = ['CFD Surrogate', 'CFD Simulation'],value = 'CFD Surrogate')
    br_options.t_final = widgets.BoundedFloatText(value = 100.0)

    return br_options


def test_run_pretreatment(build_fs_options, build_pt_options):
    fs_options = build_fs_options
    pt_options = build_pt_options

    run_pretreatment(notebook_dir, params_filename, fs_options, pt_options)

    truth_values = {'fis_0': 0.31765314961287994,
                    'conv': 0.028073229915110083,
                    'X_X': 0.2575180632106229,
                    'X_G': 0.40297527098473657,
                    'rho_x': 3.4623134078020046,
                    'rho_f': 0.0004599638814971187}

    test_values = yaml_to_dict(params_filename)

    for key, val in truth_values.items():
        assert key in test_values['pretreatment_output']
        assert val == pytest.approx(test_values['pretreatment_output'][key])


def test_run_enzymatic_hydrolysis(build_eh_options):
    eh_options = build_eh_options

    run_enzymatic_hydrolysis(notebook_dir, params_filename, eh_options, False)
    

# def test_run_bioreactor(build_br_options):
#     br_options = build_br_options
    
#     run_bioreactor(notebook_dir, params_filename, br_options, False)
