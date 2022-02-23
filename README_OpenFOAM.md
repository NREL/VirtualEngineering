# Building OpenFOAM with Submodules

When you first clone the VirtualEngineering repo, folders named `openfoam-dev` and `thirdparty-dev` are created within the `submodules` directory, but are empty by default.  To populate them and build your own version of OpenFOAM, follow these steps.

1. From the root level of the VirtualEngineering directory, do:

```bash
git submodule init # if not done already
git submodule update submodules/OpenFOAM-dev
git submodule update submodules/ThirdParty-dev
```

which will pull in the necessary files from the linked sources (a custom version of openfoam with NREL contributions along with thirdparty).

2. Load all required modules with: 

```bash
module use /nopt/nrel/ecom/hpacf/utilities/modules
module load bison/3.6.4
module load flex/2.5.39

module use /nopt/nrel/apps/modules/centos74/modulefiles
module load gcc/7.3.0 
module load openmpi/1.10.7/gcc-7.3.0
```

3. Change your directory to the OpenFOAM folder.

>**WARNING**: Note that OpenFOAM requires that the *full* path to the directory be entered in this step.  So rather than simply typing `cd submodules/OpenFOAM-dev/` as you normally would, enter the full path `cd /lustre/eaglefs/.../VirtualEngineering/submodules/OpenFOAM-dev/`.  You can check that you've done this correctly by running `pwd` and confirming the path to your current location is prefixed with `/lustre/eaglefs/`

```bash
cd /lustre/eaglefs/.../VirtualEngineering/submodules/OpenFOAM-dev/
```

4. Source the OpenFOAM bashrc file by with:

```bash
source etc/bashrc
```

5. Set the number of processors to use for compilation, e.g., 36, with:

```bash
export WM_NCOMPPROCS=36
```

6. Still within the OpenFOAM folder (and after verifying the full path is available) compile with:

```bash
./Allwmake
```

which will build all of OpenFOAM and thirdparty.  This process can take over an hour, so be sure to request an interactive session or other compute resources with a sufficient upper time limit (2 hours should suffice with `WM_NCOMPPROCS=36`).

