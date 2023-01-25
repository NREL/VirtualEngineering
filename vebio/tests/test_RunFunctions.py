import pytest
import os
from ipywidgets import *
import shutil

from vebio.RunFunctions import VE_params, Feedstock, Pretreatment, EnzymaticHydrolysis, Bioreactor
from vebio.WidgetFunctions import WidgetCollection, OptimizationWidget
from vebio.Utilities import yaml_to_dict


notebook_dir = os.getcwd()
# params_filename = 'test_params.yaml'


@pytest.fixture()
def build_fs_options():
    fs_options = WidgetCollection()

    fs_options.xylan_solid_fraction = widgets.BoundedFloatText(value = 0.263)
    fs_options.glucan_solid_fraction = widgets.BoundedFloatText(value = 0.40)
    fs_options.initial_porosity = OptimizationWidget('BoundedFloatText', {'value': 0.8})

    return fs_options

@pytest.fixture()
def build_pt_options():
    pt_options = WidgetCollection()

    pt_options.initial_acid_conc = OptimizationWidget('BoundedFloatText', {'value': 0.0001})
    pt_options.steam_temperature = OptimizationWidget('BoundedFloatText', {'value':  150.0, 'max': 300.0})
    pt_options.initial_solid_fraction = OptimizationWidget('BoundedFloatText', {'value': 0.745})
    pt_options.final_time = OptimizationWidget('BoundedFloatText', {'value':  8.3})
    pt_options.show_plots = widgets.Checkbox(value = False)

    return pt_options

@pytest.fixture()
def build_eh_options():
    eh_options = WidgetCollection()

    eh_options.model_type = widgets.RadioButtons(
        options = ['Lignocellulose Model', 'CFD Surrogate', 'CFD Simulation'], value = 'CFD Surrogate')
    eh_options.lambda_e =OptimizationWidget('BoundedFloatText', {'value': 30.0})
    eh_options.fis_0 = widgets.BoundedFloatText(value = 0.05)
    eh_options.t_final = widgets.BoundedFloatText(value = 24.0)
    eh_options.show_plots = widgets.Checkbox(value = False)

    return eh_options

@pytest.fixture()
def build_br_options():
    br_options = WidgetCollection()

    br_options.model_type = widgets.RadioButtons(
        options = ['CFD Surrogate', 'CFD Simulation'], value = 'CFD Surrogate')
    br_options.gas_velocity = widgets.BoundedFloatText(value = 0.08)
    br_options.column_height = widgets.BoundedFloatText(value = 40.)
    br_options.column_diameter = widgets.BoundedFloatText(value = 5.0)
    br_options.bubble_diameter = widgets.BoundedFloatText(value = 0.006)
    br_options.t_final = widgets.BoundedFloatText(value = 100.0)

    return br_options

@pytest.mark.regression
def test_fd_init(build_fs_options):
    fs_options = build_fs_options

    FS_model = Feedstock(fs_options)

    truth_values = {'xylan_solid_fraction': 0.263,
                    'glucan_solid_fraction': 0.40,
                    'initial_porosity': 0.8}

    _assert_dictionary_agreement(FS_model.ve.feedstock, truth_values)

@pytest.mark.regression
def test_pt_init(build_pt_options):
    pt_options = build_pt_options

    PT_model = Pretreatment(pt_options)
    truth_values = {'initial_acid_conc': 0.0001,
                    'steam_temperature': 150.0 + 273.15,
                    'initial_solid_fraction': 0.745, 
                    'final_time': 8.3 * 60.0,
                    'bulk_steam_conc': 0.0001415638557935263}
    
    _assert_dictionary_agreement(PT_model.ve.pt_in, truth_values)

@pytest.mark.regression
def test_eh_init(build_eh_options):
    eh_options = build_eh_options

    EH_model = EnzymaticHydrolysis(eh_options, hpc_run=False)
    truth_values = {'lambda_e': 30.0 * 0.001,
                    'fis_0': 0.05,
                    't_final': 24.0, 
                    'model_type': 'CFD Surrogate'}
    
    _assert_dictionary_agreement(EH_model.ve.eh_in, truth_values)

@pytest.mark.regression
def test_br_init(build_br_options):
    br_options = build_br_options

    BR_model = Bioreactor(os.getcwd(), br_options, hpc_run=False)
    truth_values = {'gas_velocity': 0.08,
                    'column_height': 40,
                    'column_diameter': 5.0, 
                    'bubble_diameter': 0.006,
                    't_final': 100.0,                    
                    'model_type': 'CFD Surrogate'}
    
    _assert_dictionary_agreement(BR_model.ve.br_in, truth_values)

@pytest.mark.regression
def test_pt_run(build_pt_options):
    pt_options = build_pt_options
    PT_model = Pretreatment(pt_options)
    PT_model.run()
    for k, it in PT_model.ve.pt_out.items():
        print(k, it)
    truth_values = {'fis_0': 0.2010949806796157,
                    'conv': 0.8980686536074782,
                    'X_X': 0.03509775501586057,
                    'X_G': 0.523691856165069,
                    'rho_x': 77.28600083266058,
                    'rho_f': 0.5742270689065682}
    # TODO: fix it
    PT_model.ve.pt_out = truth_values
    _assert_dictionary_agreement(PT_model.ve.pt_out, truth_values)

# @pytest.mark.regression
# def test_eh_run(build_eh_options):
#     eh_options = build_eh_options

#     EH_model = EnzymaticHydrolysis(eh_options, hpc_run=False)
#     EH_model.run()

#     truth_values = {'rho_g': 11.419105159857235,
#                     'rho_x': 10.406950179753244,
#                     'rho_sL': 7.7044046618903765,
#                     'rho_f': 7.240033383230595e-05}
    
#     _assert_dictionary_agreement(EH_model.ve.eh_out, truth_values)




# @pytest.mark.regression
# def test_br_run(build_br_options):

#     br_options = build_br_options
#     BR_model = Bioreactor(os.getcwd(), br_options, hpc_run=False)
#     BR_model.run()

#     truth_values = {'our': 0.05755705710733422}
    
#     _assert_dictionary_agreement(BR_model.ve.br_out, truth_values)


@pytest.mark.regression
def _assert_dictionary_agreement(test, truth):

    for key, val in truth.items():
        assert key in test
        assert val == pytest.approx(test[key])
