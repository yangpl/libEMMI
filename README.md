# libEMMI
Library for 3D controlled-source electromagnetic modelling and inversion

Author: Pengliang Yang, Harbin Institute of Technology, China

Email: ypl.2100@gmail.com

Programming language: C

Result visualization: python3, Madagascar (Scons)

System requirements: Linux OS

Compilation requirement: gcc compiler, mpicc compiler and make

Software required: fftw3

=======================================================================

1) run_adjsrc_freq2time:
--------------------------------------------

 A reproducible example to convert adjoint source from limited number of frequencies to long time series.

To run this example:

1. compile the code by running: make

2. visualize the results: python3 plot_basis_function.py



2) run_create_model: 
--------------------------------------------

A reproducible example to create resistivity models for land and marine surveys

You need to modify the parameters in main.c and then type 'make' to compile the code. Then, type './main' to run the executable to generate the geoelectric models. You can probabaly use different tools like matlab or Madagascar to visualize your model. The scripts for visualization is placed in this folder also.





3) run_fwi_land_iter30:
--------------------------------------------

A reproducible example to do full waveform inversion of 3D CSEM data with two-block model

To run this example, one has to follow the following steps:

* Go to /src and compile the source code using MPI compiler

cd /sr;

make

After compilation, the executable will be placed in /bin named as fdtd.

* Go to run_fwi_land_iter30, modify input parameters in run.sh and then launch the code under parallel computing environment.  To run it, type:

bash run.sh

By default, run.sh sets 16 sources using 16 independent MPI process.

* After running, the code will create many files about the data and the model:

------------------------------
ASCII files emf_0001.txt--emf_0016.txt store the observed data using the true resistivity model. These are created using 'mode=0' for modelling jobs.

ASCII files syn_0001.txt--syn_0016.txt store the synthetic data based on the updated model at each iteration during the inversion. After 30 iterations, they will be the synthetic data obtained from the final inverted resistivity.

ASCII files sigmisfit_0001.txt--sigmisfit_0016.txt store the significant misfit at the receiver locations in each iteration. These can be very useful to check the statistics/histogram of data error distribution.


In addition to these files, I created 3 python scripts to visualize the result. These are:

* plot_histogram.py: it will plot the histogram of data misfit.

* plot_emdata.py: it will plot the amplitude and the phase of the EM data 

* plot_scatter_sigmisfit.py: it will plot the scatter plot with color indicating the level of misfit at each receiver loation. You can change the input file in the script with other Tx index to check the scatter plot of the significant misfit associated with different sources.

The above helps you to reproduce exactly what I presented in the paper. I am using Madagascar to plot my 3D resistivity models based on SConstruct script. You would need to learn this tool by yourself if you wish to do so. The final inverted model is param_final in binary format. You can also visualize it by other tools you like.