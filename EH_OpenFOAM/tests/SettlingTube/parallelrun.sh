rm -rf 0
cp -r 0.org 0
blockMesh
transformPoints -rotate "((1 0 0)(0 0 1))"
transformPoints -scale "(0.001 0.001 0.001)"
decomposePar
srun -n 32 EHFoam -parallel
