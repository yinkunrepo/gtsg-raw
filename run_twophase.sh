#!/bin/bash

# combine two primitive surfaces and calculate the melting curve.

top=`pwd`
casedir=${top}/twophase
if [ ! -d $casedir ]; then
    mkdir -p $casedir
else
    echo "[Error] You already have an old folder, remove it first."
    exit
fi

againdir=${casedir}/liquid_again

cd $casedir
cp $top/lib/gibbs_surf.py .
cp $top/lib/2p_*.py .

cp ${top}/solid/dat_abc.txt solid_dat_abc.txt
cp ${top}/solid/dat_eos.txt solid_dat_eos.txt
cp ${top}/liquid/dat_abc.txt liquid_dat_abc.txt
cp ${top}/liquid/dat_eos.txt liquid_dat_eos.txt

python 2p_01.py

cp -r ${top}/liquid ${againdir}
cp ${casedir}/liquid_delta_S.txt ${againdir}/refS.txt
cd ${againdir}
python 1p_03.py

cd $casedir
cp ${againdir}/dat_abc.txt liquid_dat_abc.txt
cp ${againdir}/dat_eos.txt liquid_dat_eos.txt

python 2p_02.py
python 2p_03.py
python 2p_04.py > cross_points.txt
python 2p_05.py
cd $top
