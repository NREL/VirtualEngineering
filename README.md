# Virtual Engineering

This repository contains all the tools needed to launch interactive Jupyter Notebooks which control the execution of various unit operations for the end-to-end simulation of low-temperature conversion of biomass to fuel.


### Create Conda Environment

The first step is to create a Conda environment using the included `environment.yaml` file.  To create this environment, simply open a terminal, navigate to the root level of the Virtual Engineering directory and run
    
`conda env create -f environment.yaml`

which creates a Conda environment named `virteng-env` from the complete recipe of all the necessary packages enumerated in the yaml file.
    
    
### Run Notebooks

Once the Conda environment is successfully created, activate it from the command line with

`conda activate virteng-env`

then, launch the Jupyter notebook by running

`jupyter notebook`

When the notebook landing page opens, you'll be able to open the [virtual_engineering_notebook](virtual_engineering_notebook.md) to get started.
