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

1) run_adjsrc_freq2time:

 A reproducible example to convert adjoint source from limited number of frequencies to long time series.

To run this example:

1. compile the code by running: make

2. visualize the results: python3 plot_basis_function.py

2) run_create_model: 

A reproducible example to create resistivity models for land and marine surveys


3) run_fwi_land_iter30:

A reproducible example to do full waveform inversion of 3D CSEM data with two-block model

To run this example, one has to follow the following steps:

* go to /src and compile the source code using MPI compiler

cd /sr;

make

After compilation, the executable will be placed in /bin named as fdtd.

* go to run_fwi_land_iter30, modify input parameters in run.sh and then launch the code under parallel computing environment.  To run it, type:

bash run.sh

By default, run.sh sets 16 sources using 16 independent MPI process.