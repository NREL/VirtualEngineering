Getting Started
===============

Building Conda Environment
--------------------------

This Virtual Engineering repository contains all the tools needed to create Jupyter Notebooks to control the execution of various unit operations for the beginning-to-end simulation of the low-temperature conversion of lignocellulosic biomass to a fuel precursor.  The first step is to create a Conda environment using the included ``environment.yaml`` file. To create this environment, simply open a terminal, navigate to the root level of the Virtual Engineering directory and run::

	conda env create -f environment.yaml

which creates a Conda environment named ``virteng-env`` (you can change this as you see fit) using the packages enumerated in the yaml file. Additionally, it pip installs the ``vebio`` Python package which adds functionality required by the Jupyter Notebook workflow.  Once the Conda environment is successfully created, activate it from the command line with::

	conda activate virteng-env

then, launch the Notebook by running::

	jupyter notebook

When the Notebook landing page opens, you'll be able to open the ``virtual_engineering_notebook.nd`` which will allow you to launch simulations and move on to the next step.

Running Notebook
----------------

<Demo of running notebook with screenshots>
