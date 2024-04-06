#!/usr/bin/bash

echo "
mode=1 #0=modelling; 1=inversion

==================== modelling setup ================
fsrc=sources.txt
frec=receivers.txt
fsrcrec=src_rec_table.txt
frho_h=rho_init
frho_v=rho_init
chsrc=Ex
chrec=Ex,Ey
x1min=-9000
x1max=9000
x2min=-9000
x2max=9000
x3min=0
x3max=3000
n1=91
n2=91
n3=91
d1=200
d2=200
d3=25
freqs=0.25,0.75,2.25


================= grid setup ====================
nugrid=1
fx3nu=x3nu
fibathy=ibathy

================ inversion setup ===============
niter=30
nls=5
npair=5
npar=2
preco=1
bound=1
idxpar=1,2
minpar=1.0,1.0
maxpar=70.0,70.0

gamma1=2e3
gamma2=1e2

addnoise=1
offset_start=1000

">inputpar.txt

#export OMP_NUM_THREADS=4
mpirun -n 1 ../bin/fdtd $(cat inputpar.txt) 

