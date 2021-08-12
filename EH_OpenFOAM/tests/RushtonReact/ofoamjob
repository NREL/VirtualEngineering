#!/bin/bash
#SBATCH --job-name=ehreact
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --time=48:00:00
#SBATCH --account=vebio
#SBATCH --output=log.out

#required modules
#========================================
export TMPDIR=${LOCAL_SCRATCH}
module purge
module use /nopt/nrel/apps/modules/centos74/modulefiles
module load gcc/7.3.0 
module load openmpi/1.10.7/gcc-7.3.0

#source openfoam's bashrc - change accordingly
source /projects/bpms/openfoam/OpenFOAM-dev/etc/bashrc
module load conda
module load paraview

#presteps
#========================================
rm ./blockMeshDict_reactor
rm -r 0
cp -r 0.org 0
python system/write_bmesh_file.py
blockMesh -dict ./blockMeshDict_reactor
stitchMesh -perfect -overwrite inside_to_hub inside_to_hub_copy
stitchMesh -perfect -overwrite hub_to_rotor hub_to_rotor_copy
#convert to cms
transformPoints -scale "(0.01 0.01 0.01)"

#run the case
#========================================
decomposePar
srun -n 32 EHFoam -parallel
reconstructPar -newTimes
EHFoam -postProcess -func "grad(U)"
pvpython torq.py "lateralWall"
