rm ./blockMeshDict_reactor
rm -r 0
cp -r 0.org 0
python system/write_bmesh_file.py
blockMesh -dict ./blockMeshDict_reactor
stitchMesh -perfect -overwrite inside_to_hub inside_to_hub_copy
stitchMesh -perfect -overwrite hub_to_rotor hub_to_rotor_copy
#convert to cms
transformPoints -scale "(0.01 0.01 0.01)"
decomposePar
srun -n 32 EHFoam -parallel
reconstructPar -newTimes
EHFoam -postProcess -func "grad(U)"
pvpython torq.py "lateralWall"
