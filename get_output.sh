#!/usr/bin/env bash
# get output files

top=`pwd`
outdir=$top/output
if [ ! -d $outdir ]; then
    mkdir $outdir
fi

files="kechin_fit_params_with_units.txt
kechin_fit.pdf
"
for fname in $files; do
    this_file=twophase/$fname
    if [ -f $this_file ]; then
	cp $this_file $outdir
    fi
done
