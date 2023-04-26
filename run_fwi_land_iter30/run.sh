#!/usr/bin/bash

echo "
=============== modelling setup =====================
mode=1

fsrc=sources.txt
frec=receivers.txt
fsrcrec=src_rec_table.txt

frho11=rho11_init
frho22=rho22_init
frho33=rho33_init
chsrc=Ex
chrec=Ex,Ey
x1min=-3000
x1max=3000
x2min=-3000
x2max=3000
x3min=0
x3max=1500
n1=101
n2=101
n3=61
d1=60
d2=60
d3=25
freqs=0.25,1

================= inversion setup ====================
niter=30
npar=1
preco=1
bound=1
idxpar=1
minpar=1
maxpar=50.0

r1=2
r2=2
r3=2
repeat=3




addnoise=1
offset_start=600

">inputpar.txt

export OMP_NUM_THREADS=4
mpirun -n 16 ../bin/fdtd $(cat inputpar.txt) #>&out

