#!/bin/bash
#SBATCH --account=vebio
#SBATCH --time=12:00:00
#SBATCH --qos=high
#SBATCH --job-name=reactor
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=36
#SBATCH --output=log.out


module purge
module load openmpi/1.10.7/gcc-7.3.0
module load gcc
module load paraview

# Run the bioreactor simuilation
source /projects/vebio/hsitaram/VirtualEngineering/submodules/OpenFOAM-dev/etc/bashrc
. ./presteps.sh
srun -n $SLURM_NTASKS reactingTwoPhaseEulerFoam -parallel

# Run the post-processing script
source /projects/bpms/openfoam/OpenFOAM-dev/etc/bashrc
rm -rf 0.0/
rm -rf 0/
reconstructPar -latestTime
mv 0/ 0.0/
pvpython pv_extract_analyze_script.py 
