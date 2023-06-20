# GTSG

Derivation of melting curve from Gibbs thermodynamic surface for davemaoite.

## Requirement
- Bash shell
- Python3 (with	`numpy`, `scipy`, `matplotlib` packages installed)

## Run
Change to the top directory of this code and execute `./GTSG` in the terminal.

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

## Example of output files

You will find the output files in `output` folder.
- `kechin_fit_params_with_units.txt` Parameters of Kechin equation after fitting are given in this file.
- `kechin_fit.pdf` A figure shows the fitted melting curve in pressure-temperature diagram.

The `output_example` folder is just for your reference.


## Example of standard output

You will see the following text print on screen.
```shell
**************************************
Gibbs Thermodynamic Surface Generator
      _____ _______ _____  _____
     / ____|__   __/ ____|/ ____|
    | |  __   | | | (___ | |  __
    | | |_ |  | |  \___ \| | |_ |
    | |__| |  | |  ____) | |__| |
     \_____|  |_| |_____/ \_____|
           a raw version
**************************************
calculating primitive surface for solid phase...done
calculating primitive surface for liquid phase...done
calculating derived surface for two phases...done
generating output...done
cleaning temporary folders...done
```

## References

To learn more about the theory behind this code, the following papers may help.

Gibbs, J.W., 1873. *A method of geometrical representation of the thermodynamic properties of substances by means of surfaces.* Trans. Conn. Acad. 2, 382.

Yin, K., Lu, X., Zhou, H., Sun, Y., 2018. *Thermodynamic stability limit of the crystalline state from the Gibbs perspective.* Physical Review B 98, 144113. https://doi.org/10.1103/PhysRevB.98.144113
