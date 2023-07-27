#!/bin/bash
#SBATCH --account=vebio
#SBATCH --time=48:00:00
#SBATCH --job-name=sweeps
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=36
#SBATCH --output=log.out
#SBATCH --mail-user=andrew.glaws@nrel.gov
#SBATCH --mail-type=ALL


module purge
module load openmpi/1.10.7/gcc-8.4.0
#module load openmpi/1.10.7/gcc-7.3.0
module load gcc
#module load paraview
module load paraview/5.6.0

# Run the bioreactor simuilation
source /projects/vebio/aglaws/VirtualEngineering/submodules/OpenFOAM-dev/etc/bashrc
. ./presteps.sh
srun -n $SLURM_NTASKS reactingTwoPhaseEulerFoam -parallel

# Run the post-processing script
rm -rf 0.0/
rm -rf 0/
#reconstructPar -latestTime
reconstructPar -withZero -newTimes
mv 0/ 0.0/
pvpython pv_extract_analyze_script.py 
pvpython pv_extract_our.py
