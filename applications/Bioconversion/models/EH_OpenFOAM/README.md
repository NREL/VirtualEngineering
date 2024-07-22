# Running the Enzymatic Hydrolysis OpenFOAM Simulation

1. Compile OpenFOAM-dev using the instructions found in [submodules/README.md](../submodules/README.md).  If OpenFOAM has already been installed using these instructions, skip to Step 2. 

2. Load the modules necessary to make the enzymatic hydrolysis examples:

```bash
module use /nopt/nrel/apps/modules/centos74/modulefiles
module load gcc/7.3.0 
module load openmpi/1.10.7/gcc-7.3.0
```

3. From the root level of the VirtualEngineering directory, source the OpenFOAM bashrc file with:

```bash
source submodules/OpenFOAM-dev/etc/bashrc
```

substituting the correct path to the bashrc file as needed if the location of your OpenFOAM installation differs.

4. Change directory to the `EHFoam` folder and do `wmake`

```bash
cd EH_OpenFOAM/EHFoam/
wmake
```

5. Try out any of the cases in the `tests` directory, for example, `RushtonReact` is a production case for lignocellulose digestion. 

* The biomass inputs are in constant/globalVars. 
* constant/EHProperties has transport coefficients and some solver knobs/parameters. 
* constant/MRFproperties has the rotational speed

6. The coupled runs creates a file called `integrated_quantities.dat` which 
has integrated values of solids volume fractions and averaged dissolved sugar concentrations.
This file can be postprocessed to obtain conversion data
