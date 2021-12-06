Using OpenFOAM Unit Models
==========================

Building OpenFOAM
-----------------

When you first clone the VirtualEngineering repo, folders named ``openfoam-dev`` and ``thirdparty-dev`` are created within the submodules directory, but are empty by default. To populate them and build your own version of OpenFOAM, follow these steps.  From the root level of the VirtualEngineering directory, do::

    git submodule init
    git submodule update

which will pull in the necessary files from the linked sources (a custom version of openfoam with NREL contributions along with thirdparty).  Next, load all required modules with::

    module use /nopt/nrel/ecom/hpacf/utilities/modules
    module load bison/3.6.4
    module load flex/2.5.39

    module use /nopt/nrel/apps/modules/centos74/modulefiles
    module load gcc/7.3.0 
    module load openmpi/1.10.7/gcc-7.3.0

Once these modules are loaded, change your directory to the OpenFOAM folder::

    cd /lustre/eaglefs/.../VirtualEngineering/submodules/OpenFOAM-dev/

.. warning:: Note that OpenFOAM requires that the full path to the directory be entered in this step. So rather than simply typing cd submodules/OpenFOAM-dev/ as you normally would, enter the full path cd /lustre/eaglefs/.../VirtualEngineering/submodules/OpenFOAM-dev/. You can check that you've done this correctly by running pwd and confirming the path to your current location is prefixed with /lustre/eaglefs/.

Source the OpenFOAM bashrc file by with::

    source etc/bashrc

and set the number of processors to use for compilation, e.g., 36, with::

    export WM_NCOMPPROCS=36

Still within the OpenFOAM folder (and after verifying the full path is available) compile with::

    ./Allwmake

which will build all of OpenFOAM and thirdparty. This process can take over an hour, so be sure to request an interactive session or other compute resources with a sufficient upper time limit (2 hours should suffice with ``WM_NCOMPPROCS=36``).

Running OpenFOAM Simulations
----------------------------

To run the Enzymatic Hydrolysis OpenFOAM Simulation, first compile OpenFOAM-dev using the instructions found in :ref:`Building OpenFOAM`. If OpenFOAM has already been installed using these instructions, proceed with loading the modules necessary to make the enzymatic hydrolysis examples::

    module use /nopt/nrel/apps/modules/centos74/modulefiles
    module load gcc/7.3.0 
    module load openmpi/1.10.7/gcc-7.3.0

From the root level of the VirtualEngineering directory, source the OpenFOAM bashrc file with::

    source submodules/OpenFOAM-dev/etc/bashrc

substituting the correct path to the bashrc file as needed if the location of your OpenFOAM installation differs.  Next, change directory to the EHFoam folder and do ``wmake``::

    cd EH_OpenFOAM/EHFoam/
    wmake

Try out any of the cases in the tests directory, for example, RushtonReact is a production case for lignocellulose digestion. A few orientation notes on the specification files:

* The biomass inputs are in constant/globalVars.
* constant/EHProperties has transport coefficients and some solver knobs/parameters.
* constant/MRFproperties has the rotational speed

The coupled runs creates a file called integrated_quantities.dat which has integrated values of solids volume fractions and averaged dissolved sugar concentrations. This file can be postprocessed to obtain conversion data
