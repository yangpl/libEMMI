#include "cstd.h"

#ifndef PI
#define PI (3.1415926535897932384)
#endif

void adjsrc_freq2time_init(float omega0, float dt, float *omegas, int nfreq, int nt, int ntr,
			   float tol_, float eps_, int nrestart_, int niter_, int verb_);
void adjsrc_freq2time_close();
void adjsrc_init(float *mm_, float *dd_);
void adjsrc_freq2time();

int main(void)
{
  int ifreq,i,j,ntrace;
  float dt=0.004;
  int nt=1000;
  float f0=1.;
  float omega0 = 2.*PI*f0;
  const int nfreq = 3;
  float freqs[3] = {0.25, 0.75, 1.25};
  float omegas[3];
  float *xx, *yy;

  float adjsrc_tol = 1e-6;
  float adjsrc_eps = 1e-3;
  int adjsrc_nrestart = 30;
  int adjsrc_niter = 100;
  int verb = 1;

  for(ifreq=0; ifreq<nfreq; ++ifreq) omegas[ifreq] = 2.*PI*freqs[ifreq];

  ntrace = 2*nfreq;
  xx = alloc1float(nt*ntrace);
  yy = alloc1float(2*nfreq*ntrace);
  memset(xx, 0, nt*ntrace*sizeof(float));
  for(j=0; j<ntrace; j++){
    for(ifreq=0; ifreq<2*nfreq; ++ifreq) {
      i = ifreq + 2*nfreq*j;
      yy[i] = (ifreq==j)?1.:0.;
    }
  }

  adjsrc_freq2time_init(omega0, dt, omegas, nfreq, nt, ntrace, adjsrc_tol, adjsrc_eps, adjsrc_nrestart, adjsrc_niter, verb);
  adjsrc_init(xx, yy);
  adjsrc_freq2time();/* most computationally intensive part hidden here */    
  adjsrc_freq2time_close();      

  FILE *fp = fopen("basis_function", "wb");
  fwrite(xx, nt*ntrace*sizeof(float), 1, fp);
  fclose(fp);


  free(xx);
  free(yy);

}
