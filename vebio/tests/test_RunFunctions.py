import pytest
import os
from ipywidgets import *
import shutil

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


@pytest.mark.regression
def test_run_pt(build_fs_options, build_pt_options):
    fs_options = build_fs_options
    pt_options = build_pt_options

    run_pretreatment(notebook_dir, params_filename, fs_options, pt_options)

    test_values = yaml_to_dict(params_filename)

    truth_values = {'fis_0': 0.31765314961287994,
                    'conv': 0.028073229915110083,
                    'X_X': 0.2575180632106229,
                    'X_G': 0.40297527098473657,
                    'rho_x': 3.4623134078020046,
                    'rho_f': 0.0004599638814971187}

    _assert_dictionary_agreement(test_values['pretreatment_output'], truth_values)

@pytest.mark.regression
def test_run_eh_surrogate(build_eh_options):
    eh_options = build_eh_options

    run_enzymatic_hydrolysis(notebook_dir, params_filename, eh_options, False)

    test_values = yaml_to_dict(params_filename)

    truth_values = {'rho_g': 11.419105159857235,
                    'rho_x': 10.406950179753244,
                    'rho_sL': 7.7044046618903765,
                    'rho_f': 7.240033383230595e-05}

    _assert_dictionary_agreement(test_values['enzymatic_output'], truth_values)

@pytest.mark.regression
def test_run_br_surrogate(build_br_options):
    br_options = build_br_options
    
    run_bioreactor(notebook_dir, params_filename, br_options, False)

    test_values = yaml_to_dict(params_filename)

    truth_values = {'our': 0.05755705710733422}

    _assert_dictionary_agreement(test_values['bioreactor_output'], truth_values)

@pytest.mark.regression
def _assert_dictionary_agreement(test, truth):

    for key, val in truth.items():
        assert key in test
        assert val == pytest.approx(test[key])
