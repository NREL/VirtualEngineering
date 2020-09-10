#!/bin/bash
NPROC=$1
PROBLEMNAME=paddle
echo ${PROBLEMNAME} > input_genmap
echo 0.01 >> input_genmap
export NEK_SOURCE_ROOT=../submodules/Nek5000
export PATH=${NEK_SOURCE_ROOT}/bin:${PATH}
genmap < input_genmap
echo ${PROBLEMNAME} > SESSION.NAME
echo $PWD >> SESSION.NAME
makenek paddle -y
mpirun -n ${NPROC} ./nek5000
