# <center> MolPG: detect point group of a molecule </center>

Usage: `MolPG xxx.gjf/xxx.xyz`

Supports Gaussian `.gjf` format or standard `.xyz` format.
There is a tolerance value `tol` in `main.cpp` you may change.

To compile yourself, you neead to install **Eigen** 3.4+, and **fmt** 8.0+ first.
A C++ compiler with C++17 standard supported is required.

Some ideas come from reference [10.1002/jcc.23493](https://onlinelibrary.wiley.com/doi/10.1002/jcc.23493).

