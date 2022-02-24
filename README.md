# Virtual Engineering

The Virtual Engineering (VE) project provides a set of Python modules to enable users to connect mathematical models of unit operations to predict outcomes for an entire chemical process. This VE approach is currently being developed to support the beginning-to-end simulation of the low-temperature conversion of lignocellulosic biomass to a fuel precursor. More details can be found in the [VE documentation](https://virtualengineering.readthedocs.io/en/latest/index.html).

![Test Status](https://github.com/NREL/VirtualEngineering/actions/workflows/test_vebio.yml/badge.svg)
[![Documentation Status](https://readthedocs.org/projects/virtualengineering/badge/?version=latest)](https://virtualengineering.readthedocs.io/en/latest/?badge=latest)


## Getting Started

New users are encouraged to review the [Getting Started](https://virtualengineering.readthedocs.io/en/latest/how_to_guides/getting_started.html#getting-started) guide which describes how to create the Conda environment and run the Jupyter Notebook.

## Developer Quick Start

First, create a Conda environment from the `environment.yaml` file by opening a terminal, navigating to the root level of the Virtual Engineering directory, and running

`conda env create -f environment.yaml -n <env_name>`

where `<env_name>` is the desired name of your VE Conda environment. This Conda installs the required packages and pip installs the `vebio` Python package.  Next, activate the environment with

`conda activate <env_name>`

and launch the Notebook by running

`jupyter notebook`

When the Notebook landing page opens, you'll be able to open the [virtual_engineering_notebook](virtual_engineering_notebook.md) which will guide you through the next steps.

## Submodules

Git has a mechanism for code in a repository to depend on code in another repository. The dependency is known as a [submodule](https://git-scm.com/book/en/v2/Git-Tools-Submodules). Submodules are specified in `.gitmodules` (a hidden file on some systems) and are located in subfolders of `submodules/` of this repository. Initiate the submodule system with

`git submodule init`

Then update specific modules with 

`git submodule update submodules/[name of submodule]`

This command will perform a git clone or pull for that module. To run the enzymatic hydrolysis CFD model, update `Nek5000`. To run the well-mixed "Ligncellulose" model (no CFD), update`CEH_EmpiricalModel`. To run bioreactor CFD simulation, update `OpenFOAM-dev` and `ThirdParty-dev`. Also see [README_OpenFOAM.md](README_OpenFOAM.md).

