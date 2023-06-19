# GTSG

Derivation of melting curve from Gibbs thermodynamic surface for davemaoite.

## Requirement
- Bash shell
- Python3 (with	`numpy`, `scipy`, `matplotlib` packages installed)

## Run
Execute `./GTSG` in the terminal.

## Description of file structure

### `input` folder

- `solid_input.txt` file
- `liquid_input.txt` file

The above two files provide data of volume (given by lattice ratio
with respect to a given basis), pressure, energy and temperature.
In this example (i.e., davemaoite), these data are computed
by VASP software based on a large number of NVT-ensemble
first principles molecular dynamics (FPMD) simulations.

### `lib` folder

This folder contains all python codes. `1p_xx.py` are codes for
calculating one-phase surface. `2p_xx.py` are codes for calculating
melting curve.

### `output_example` folder
After run, you will find the `output` folder. This `output_example` folder is just
for your reference.

- `kechin_fit_params_with_units.txt` Parameters of Kechin equation after fitting are given in this file.
- `kechin_fit.pdf` A figure shows the fitted melting curve in pressure-temperature diagram.
