# Virtual Engineering Notebook

This repository contains all the demos and development tools needed to launch interactive Jupyter notebooks which control the execution of external unit operations.  The first step is to create a Conda environment using one of the following two methods.

1. **Environment File (Recommended)**
    To create a Conda environment from the environment.yaml file, simply open a terminal and run
    
    `conda env create -f environment.yaml`
    
    which creates a Conda environment named `virteng` from the complete list of all the necessary virtual engineering packages enumerated in the yaml file. 

2. **Command Line**

    For a more customizable option, the Conda environment can be created by specifying the desired packages directly from the command line

    `conda create -n virteng python=3.7 jupyter pyyaml numpy matplotlib scipy xlrd`

    which creates a conda environment named `virteng` containing the listed packages.  Feel free to add, remove, or modify these packages to suit your specific virtual engineering needs.

Once the Conda environment is successfully created, activate it with

`conda activate virteng`

then, launch the jupyter notebook with

`jupyter notebook`

When the notebook landing page opens, you'll be able to open the [VirtualEngineeringDemo.ipynb](demo/VirtualEngineeringDemo.ipynb) to get started.
