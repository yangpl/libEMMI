# libEMMI
library for 3D controlled-source electromagnetic modelling and inversion

Author: Pengliang Yang, Harbin Institute of Technology, China

Email: ypl.2100@gmail.com

Programming language: C

Result visualization: python3, Madagascar (Scons)

System requirements: Linux OS

Compilation requirement: gcc compiler, mpicc compiler and make

Software required: fftw3
=======================================================================

* run_adjsrc_freq2time:
 A reproducible example to convert adjoint source from limited number of frequencies to long time series.

To run this example:

1. compile the code by running: make

2. visualize the results: python3 plot_basis_function.py

* run_create_model: 
A reproducible example to create resistivity models for land and marine surveys


* run_fwi_land:
A reproducible example to do full waveform inversion of 3D CSEM data with two-block model

