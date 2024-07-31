# libEMMI
Library for 3D controlled-source electromagnetic modelling and inversion

Author: Pengliang Yang, Harbin Institute of Technology, China

Email: ypl.2100@gmail.com

Programming language: C, Fortran

Result visualization: python3, Madagascar (Scons) and gnuplot

System requirements: Linux OS

Compilation requirement: gcc compiler, mpicc compiler and make

Software required: fftw3 (https://fftw.org/)

## Credits:

1. Pengliang Yang, 3D fictitious wave domain CSEM inversion by adjoint source estimation, 2023 Computers & Geosciences , Vol. 180 p. 105441
[doi:10.1016/j.cageo.2023.105441](https://doi.org/10.1016/j.cageo.2023.105441)

2. Pengliang Yang, libEMM: A fictious wave domain 3D CSEM modelling library bridging sequential and parallel GPU implementation, 2023 Computer Physics Communications, Vol. 288 p. 108745 [doi:10.1016/j.cpc.2023.108745](https://doi.org/10.1016/j.cpc.2023.108745)

3. Pengliang Yang and Rune Mittet, Controlled-source electromagnetics modelling using high order finite-difference time-domain method on a nonuniform grid 2023 Geophysics , Vol. 88, No. 2 Society of Exploration Geophysicists p. E53-E67
[doi:10.1190/geo2022-0134.1](https://doi.org/10.1190/geo2022-0134.1)


## Software Structure

* src: source code 
* include: header files
* bin: the directory where executable will be stored
* run_adjsrc_freq2tiime: adjoint source estimation example
* run_create_model: template to create resistivity model
* run_fwi_land_iter30: template for 3D CSEM inversion
* run_fwi_land_iter1: template to run only 1st iteration


## Examples

1) run_adjsrc_freq2time:


 A reproducible example to convert adjoint source from limited number of frequencies to long time series.

To run this example:

1. compile the code by running: make

2. visualize the results: python3 plot_basis_function.py


--------------------------------------------
2) run_create_model: 


A reproducible example to create resistivity models for land and marine surveys

You need to modify the parameters in main.c and then type 'make' to compile the code. Then, type './main' to run the executable to generate the geoelectric models. You can probabaly use different tools like matlab or Madagascar to visualize your model. The scripts for visualization is placed in this folder also.




--------------------------------------------
3) run_fwi_land_iter30:


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
The input before inversion:

You would first need to create input file. There is a Fortran routine create_acquisition.f90 requiring gfortran to compile and output the acquisition files in ASCII format. The acquisition files include:

* sources.txt, describes the source coordinates by x, y, z, azimuth, dip and iTx;

* receivers.txt, describes the receiver coordinates by x, y, z, azimuth, dip and iRx;

* src_rec_table.txt, gives the connection table showing which source uses which receivers to record EM data.

You can visualize the survey configuration by running gnuplot script:

* gnuplot plot_survey_layout.gnu

------------------------------
The output after inversion:

* ASCII files emf_0001.txt--emf_0016.txt store the observed data using the true resistivity model. These are created using 'mode=0' for modelling jobs.

* ASCII files syn_0001.txt--syn_0016.txt store the synthetic data based on the updated model at each iteration during the inversion. After 30 iterations, they will be the synthetic data obtained from the final inverted resistivity.

* ASCII files sigmisfit_0001.txt--sigmisfit_0016.txt store the significant misfit at the receiver locations in each iteration. These can be very useful to check the statistics/histogram of data error distribution.

* ASCII file iterate.txt gives all relevant information about nonlinear l-BFGS optimization. It gives the details such as:

 l-BFGS memory length: 5

 Maximum number of iterations: 30

 Convergence tolerance: 1.00e-08

 maximum number of line search: 10

 initial step length: alpha=1

 iter    fk       fk/f0      ||gk||    alpha    nls   ngrad

  0   1.65e+02  1.00e+00   1.09e+01  1.00e+00    0     0

  1   1.30e+02  7.85e-01   7.81e+00  1.00e+00    0     1

  2   1.14e+02  6.92e-01   2.24e+00  1.00e+00    0     2

  3   1.10e+02  6.63e-01   2.06e+00  1.00e+00    0     3

  4   8.24e+01  4.99e-01   1.88e+00  1.00e+00    0     4

  ...

I output all relevant information during nonlinear CSEM inversion into ASCII file out.txt by 'mpirun -n 16 ../bin/fdtd $(cat inputpar.txt) >&out.txt' in run.sh. You can check out all these details.


## Visualization of the numerical result

I created 3 python scripts to visualize the result. These are:

* plot_histogram.py: it will plot the histogram of data misfit.

* plot_emdata.py: it will plot the amplitude and the phase of the EM data 

* plot_scatter_sigmisfit.py: it will plot the scatter plot with color indicating the level of misfit at each receiver loation. You can change the input file in the script with other Tx index to check the scatter plot of the significant misfit associated with different sources.

The above helps you to reproduce exactly what I presented in the paper. I am using Madagascar to plot my 3D resistivity models based on SConstruct script. You would need to learn this tool by yourself (www.ahay.org) if you wish to do so. Running 'scons' will create a folder Fig including figures in vpl format and stored.

The final inverted model is param_final in binary format. You can also visualize it by other tools you like.

