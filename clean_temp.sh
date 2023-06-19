#!/usr/bin/env bash
# clean temporary data

for dir in solid liquid twophase; do
    if [ -d $dir ]; then
	rm -rf $dir
    fi
done
