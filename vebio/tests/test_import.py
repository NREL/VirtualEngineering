import pytest
import sys
import os

def test_pt_import_simple():
    import pt

def test_pt_import():
    sys.path.append(os.path.join(os.getcwd(), 'pretreatment_model/test'))
    os.chdir(os.path.join(os.getcwd(), 'pretreatment_model/test'))
    import pt
