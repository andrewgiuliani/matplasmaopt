# matplasmaopt
Release accompanying the submission of "Single-stage gradient-based stellarator coil design: Optimization for near-axis quasi-symmetry"

Run driver.m to run the optimization described in Problem 1 of the paper.  For portability reasons, this version does not call the optimized C++ implemenation of the Biot-Savart law, so reducing the number of quadrature points or the number of coil Fourier coefficients may be necessary for speed.
