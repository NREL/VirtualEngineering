import pytest
import os
from ipywidgets import *

from virteng.WidgetFunctions import WidgetCollection, ValueRangeWidget

@pytest.mark.unit
def test_widget_collection():
    test_options = WidgetCollection()

    assert isinstance(test_options, WidgetCollection)

    test_options.a = widgets.BoundedFloatText(value = 30.0, min = 0.0, max = 100.0)
    test_options.b = widgets.BoundedFloatText(value = 40.0, min = 0.0, max = 100.0)
    test_options.b.scaling_fn = lambda val : val + 1.0
    test_options.c = widgets.BoundedFloatText(value = 50.0, min = 0.0, max = 100.0)
    test_options.c.scaling_fn = lambda val : val * 0.01

    for widget_name, widget in test_options.__dict__.items():
        assert widget_name in ['a', 'b', 'c']
        assert hasattr(widget, 'value')
        assert hasattr(widget, 'min')
        assert hasattr(widget, 'max')

    test_dict = test_options.export_widgets_to_dict()

    truth_dict = {'a': 30.0, 'b': 41.0, 'c': 0.5}

    for key, val in truth_dict.items():
        assert key in test_dict
        assert test_dict[key] == pytest.approx(val)


@pytest.mark.unit
def test_value_range_widget():
    vrw = ValueRangeWidget('Porosity', 'The value of porosity.', [0.0, 1.0], [0.1, 0.9], 0.01)

    assert hasattr(vrw, 'lower')
    assert hasattr(vrw, 'upper')

    vrw.lower.value = -100.0

    assert vrw.lower.value == pytest.approx(0.0)

    vrw.lower.value = 100.0

    assert vrw.lower.value == pytest.approx(0.9)
    assert vrw.upper.value == pytest.approx(1.0)

    vrw.upper.value = 0.8

    assert vrw.lower.value == pytest.approx(0.8)
    assert vrw.upper.value == pytest.approx(0.9)

    vrw.upper.value = -100.0

    assert vrw.lower.value == pytest.approx(0.0)
    assert vrw.upper.value == pytest.approx(0.8)

