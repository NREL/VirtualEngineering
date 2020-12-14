# About openfoam submodules

When you clone the repo, folders called openfoam-dev and thirdparty-dev are created in submodules

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

step 3: set number of compilation procs, lets say I want to use 6 processors, then do, export WM_NCOMPPROCS=6

step 4: invoke allwmake as $./Allwmake - this will build thirdparty and all of openfoam
