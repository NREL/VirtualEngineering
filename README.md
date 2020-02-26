# Virtual Engineering

This repository contains all the demos and development tools needed to launch interactive Jupyter notebooks which control the execution of external unit operations.



### Create Conda Environment

The first step is to create a Conda environment using one of the following two methods.

1. **Environment File (Recommended)**
  
    To create a Conda environment from the environment.yaml file, simply open a terminal and run
    
    `conda env create -f environment.yaml`    

    which creates a Conda environment named `virteng` from the complete list of all the necessary virtual engineering packages enumerated in the yaml file.
    
    
    
2. **Command Line**

    For a more customizable option, the Conda environment can be created by specifying the desired packages directly from the command line

    `conda create -n virteng python=3.7 jupyter pyyaml numpy matplotlib scipy xlrd`

    which creates a conda environment named `virteng` containing the listed packages.  Feel free to add, remove, or modify these packages to suit your specific virtual engineering needs.



### Run Notebooks

Once the Conda environment is successfully created, activate it from the command line with

`conda activate virteng`

then, launch the Jupyter notebook by running

`jupyter notebook`

When the notebook landing page opens, you'll be able to open the [VirtualEngineeringDemo.ipynb](demo/VirtualEngineeringDemo.ipynb) to get started.


### about openfoam submodules

when you clone the repo, folders called openfoam-dev and thirdparty-dev are created in submodules

do

git submodule init
git submodule update

This will check out the custom version of openfoam with NREL contributions along with thirdparty.

To install openfoam on eagle you will need to first load some modules

module use /nopt/nrel/ecom/hpacf/compilers/modules
module use /nopt/nrel/ecom/hpacf/utilities/modules
module use /nopt/nrel/ecom/hpacf/binaries/modules
module use /nopt/nrel/ecom/hpacf/software/modules/gcc-7.4.0

module load bison/3.0.5
module load flex/2.5.39
module load gcc/7.3.0
module load openmpi/1.10.7/gcc-7.3.0

step 1: source openfoam bashrc as $source submodules/OpenFOAM-dev/etc/bashrc
step 2: change directory to openfoam folder : 
somehow openfoam needs the complete path to the directory, if ever you are installing on scratch or projects on eagle
do change folder to cd /lustre/eaglefs/projects/vebio/VirtualEngineering/submodules/OpenFOAM-dev and 
NOT /projects/vebio/VirtualEngineering/submodules/OpenFOAM-dev
step 3: set number of compilation procs, lets say I want to use 6 processors, then do, export WM_NUMCOMPPROCS=6
step 4: invoke allwmake as $./Allwmake - this will build thirdparty and all of openfoam

To run the bioreactor test-case
==============================

see folder VirtualEngineering/bioreactor/bubble_column
To run on eagle, you have just submit the ofoamjob script.
It takes about 6 hours to run to get to steady-state.

You can use the pv_extract_analyze_script.py script (also invoked through pvjob script) to 
analyse oxygen concentration and gas hold-up. This script uses paraview python, which is part of 
paraview module on eagle. This script writes a file called volume_avg.dat which stores the 
time history (first column) of oxygen concentration in mol/m3 (second column). 
The fifth column stores gas hold-up.

