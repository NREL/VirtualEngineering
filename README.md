# Virtual Engineering

The Virtual Engineering (VE) project provides a set of Python modules to enable users to connect mathematical models of unit operations to predict outcomes for an entire chemical process. This VE approach is currently being developed to support the beginning-to-end simulation of the low-temperature conversion of lignocellulosic biomass to a fuel precursor. More details can be found in the [VE documentation](https://virtualengineering.readthedocs.io/en/latest/index.html).

[![Documentation Status](https://readthedocs.org/projects/virtualengineering/badge/?version=latest)](https://virtualengineering.readthedocs.io/en/latest/?badge=latest)
![Test Status](https://github.com/NREL/VirtualEngineering/actions/workflows/test_vebio.yml/badge.svg)


### Getting Started

New users are encouraged to review the [Getting Started](https://virtualengineering.readthedocs.io/en/latest/how_to_guides/getting_started.html#getting-started) guide which describes how to create the Conda environment and run the Jupyter Notebook.


### Developer Quick Start

Create a Conda environment from the `environment.yaml` file. To create this environment, open a terminal, navigate to the root level of the Virtual Engineering directory, and run

`conda env create -f environment.yaml -n <env_name>`

where `<env_name>` is the desired name of the Conda environment. This Conda installs the necessary packages and pip installs the `vebio` Python package which adds functionality required by the Jupyter Notebook workflow.  Next, activate the environment with

`conda activate <env_name>`

and launch the Notebook with

`jupyter notebook`

When the Notebook landing page opens, you'll be able to open the [virtual_engineering_notebook](virtual_engineering_notebook.md) which will guide you through the next steps.
