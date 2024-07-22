# import sys
# import os
# import contextlib
# import subprocess
# import numpy as np

from virteng.Utilities import dict_to_yaml, yaml_to_dict, print_dict

# root_path = os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir))
# cwd = os.getcwd()

class VE_params(object):
    ''' This  class is used for storing Virtual Engineering parameters 
        so they can be accesed from any model. It uses the Borg pattern. 
        The Borg pattern (also known as the Monostate pattern) is a way to
        implement singleton behavior, but instead of having only one instance
        of a class, there are multiple instances that share the same state. In
        other words, the focus is on sharing state instead of sharing instance.
    '''

    __shared_state = {}

    def __init__(self):
        self.__dict__ = self.__shared_state

    @classmethod
    def load_from_file(cls, yaml_filename, verbose=False):
        ve = cls()
        for k, item in yaml_to_dict(yaml_filename, verbose).items():
            setattr(ve, k, item)
        return ve

    def write_to_file(self, yaml_filename, merge_with_existing=False, verbose=False):
        dict_to_yaml(self.__dict__,  yaml_filename, merge_with_existing, verbose)
    
    def __str__(self):
        return str(print_dict(self.__dict__))
    