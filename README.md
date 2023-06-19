# GTSG

## Description
Derivation of melting curve from Gibbs thermodynamic surface for davemaoite.

## Run
Execute `./GTSG` in the terminal.

## Requirement
Python3 (with numpy, scipy, matplotlib packages installed)

## File structure

### `input` folder
This folder contains text files as input data:
..* `solid_input.txt`
..* `liquid_input.txt`

The above two files inclue data of volume (given by lattice ratio
with respect to a given basis), pressure, energy and temperature
which are computed from NVT-ensemble FPMD simulations for solid phase
and liquid phase, respectively.

### `lib` folder
This folder contains python codes.

### shell codes
..* `run_onephase.sh`
..* `run_twophase.sh`
..* ...

### `output` folder
After run, you will find the output files in this folder.

..* `kechin_fit_params_with_units.txt` Parameters of Kechin equation after fitting are given in this file.
..* `kechin_fit.pdf` A figure shows the fitted melting curve in pressure-temperature diagram.
