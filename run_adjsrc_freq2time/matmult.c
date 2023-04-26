/* linear operator (lop) for matrix vector product
 * 
 * Note the matrix B is given on the fly based on complex frequency due to 
 * fictious domain formulation in Mittet (2010).
 *
 *   Copyright (c) 2020, Harbin Institute of Technology, China
 *   Author: Pengliang Yang
 *   E-mail: ypl.2100@gmail.com
 *   Homepage: https://yangpl.wordpress.com
 */
#include "cstd.h"

#ifdef _OPENMP
#include <omp.h>
#endif

float dt, omega0;
int nfreq, nt, ntr;
float *omegas;

void sf_matmult_init(float omega0_, float dt_, float *omegas_,int nfreq_,int nt_, int ntr_)
{
  omega0 = omega0_;
  dt = dt_;
  nfreq = nfreq_;
  nt = nt_;
  omegas = omegas_;
  ntr = ntr_;
}

void sf_adjnull (bool adj /* adjoint flag */, 
		 bool add /* addition flag */, 
		 int nx   /* size of x */, 
		 int ny   /* size of y */, 
		 float * x, 
		 float * y) 
/*< adjnull version for complex data. >*/
{
  int i;
    
  if(add) return;
    
  if(adj) {
    for (i = 0; i < nx; i++) x[i] = 0.0;
  } else {
    for (i = 0; i < ny; i++) y[i] = 0.0;
  }
}

void sf_matmult_lop (bool adj, bool add, int nx, int ny, float *x, float *y)
/*< operator >*/
{
  int it, ifreq, i;
  float bb, ss, cc;

  if(nx!=nt*ntr || ny!=2*nfreq*ntr) 
    fprintf(stderr, "input dimension error: nt!=%d or nfreq!=%d\n",nx,ny);
  sf_adjnull (adj,add,nx,ny,x,y);

#ifdef _OPENMP
#pragma omp parallel for default(none)		\
  schedule(static)				\
  private(i,it,ifreq,bb,cc,ss)			\
  shared(x,y,omegas,ntr,nt,nfreq,omega0,dt,adj)
#endif
  for(i=0; i<ntr; i++){
    for (it=0; it < nt; it++) {
      for (ifreq=0; ifreq < nfreq; ifreq++) {
	/*ny=nfreq; nx=nt; */
	bb = sqrt(omegas[ifreq]*omega0)*it*dt;
	cc = cos(bb);
	ss = sin(bb);
	bb = exp(-bb);
	if (adj) {
	  x[it + nt*i] += bb*cc*y[ifreq +       2*nfreq*i];
	  x[it + nt*i] += bb*ss*y[ifreq+nfreq + 2*nfreq*i];
	} else {
	  y[ifreq      +2*nfreq*i] += bb*cc*x[it+nt*i];
	  y[ifreq+nfreq+2*nfreq*i] += bb*ss*x[it+nt*i];
	}
      }
    }
  }

  
}
