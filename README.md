# Virtual Engineering Notebook

This repository contains all the demos and development tools needed to launch interactive Jupyter notebooks which control the execution of external pieces of code.  The first step is to create a Conda environment using the following command

`conda create -n virteng python=3.7 jupyter pyyaml numpy matplotlib`

which creates a conda environment named `virteng` with the five specified packages.  Once the Conda environment is successfully created, activate it with

`conda activate virteng`

then, launch the jupyter notebook with

`jupyter notebook`

When the notebook landing page opens, you'll be able to open the [VirtualEngineeringDemo.ipynb](VirtualEngineeringDemo.ipynb) to get started.
