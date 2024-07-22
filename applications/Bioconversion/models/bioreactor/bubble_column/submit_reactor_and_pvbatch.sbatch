#!/bin/bash
#SBATCH --account=vebio
#SBATCH --time=04:00:00
#SBATCH --partition=short
##SBATCH --qos=high
#SBATCH --job-name=reactor
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=36
#SBATCH --output=log.out


module purge
module load conda
conda activate virteng-env
python br_preprocess_script.py
conda deactivate

module purge
module load openmpi/4.1.0/gcc-8.4.0 
module load gcc
source /projects/vebio/eyoung/VirtualEngineering/submodules/OpenFOAM-dev/etc/bashrc
. ./presteps.sh
srun -n $SLURM_NTASKS reactingTwoPhaseEulerFoam -parallel

# Run the post-processing script
module load paraview/5.6.0
rm -rf 0.0/
rm -rf 0/
reconstructPar -latestTime
mv 0/ 0.0/
pvpython pv_extract_analyze_script.py 


module purge
module load conda
conda activate virteng-env
module load paraview/5.6.0
python br_postprocess_script.py