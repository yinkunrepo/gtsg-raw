#!/usr/bin/env bash

# construct the primitive surface for a single phase.

phase=$1
top=`pwd`
if [ ! -d ${phase} ]; then
    mkdir -p ${phase}
else
    echo "[Error] You already have an old folder, remove it first."
    exit
fi
cd ${phase}
cp ${top}/input/${phase}_input.txt input.txt
cp ${top}/lib/1p_*.py .
cp ${top}/lib/gibbs_surf.py .
python 1p_01.py
python 1p_02.py
echo 0 > ref_index.txt
echo 0. > refS.txt
cp mid_T.txt refT.txt
python 1p_03.py
cd $top
