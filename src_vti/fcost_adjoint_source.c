/* compute data misfit and adjoint source in frequency domain for selected offset;
 * convert adjoint source from freq to time domain via least-squares minimization. 
 *
 * Copyright (c) 2020 Pengliang Yang. All rights reserved.
 * Email: ypl.2100@gmail.com
 */
#include "cstd.h"
#include "emf.h"
#include "acqui.h"
#include "constants.h"
#include "mpi_info.h"


void write_misfit(acqui_t *acqui, emf_t *emf, char *fname, float ***obj);

void adjsrc_freq2time_init(float omega0, float dt, float *freqs, int nfreq, int nt, int ntrace,
			   float tol_, float eps_, int nrestart_, int niter_, int verb_);
void adjsrc_freq2time_close();
void adjsrc_init(float *mm_, float *dd_);
void adjsrc_freq2time();

/*----------------------------------------------*/
float adjsrc_tol, adjsrc_eps;
int adjsrc_niter, adjsrc_nrestart;


/*----------------------------------------------*/
void fcost_adjsrc_init()
{
  if(!getparint("adjsrc_niter", &adjsrc_niter)) adjsrc_niter = 100;
  if(!getparint("adjsrc_nrestart", &adjsrc_nrestart)) adjsrc_nrestart = 10;
  if(!getparfloat("adjsrc_tol", &adjsrc_tol)) adjsrc_tol = 1e-7;
  if(!getparfloat("adjsrc_eps", &adjsrc_eps)) adjsrc_eps = 1e-4;//regularization parameter
}


/*----------------------------------------------*/
float fcost_adjoint_source(acqui_t *acqui, emf_t *emf)
/*< compute cost function and the adjoint source for FWI >*/
{
  float fcost;
  int ichrec, irec, ifreq, it, ntrace, i, j;
  float _Complex *src_fd;
  float *xx, *yy;
  float s;

  /* 0. compute basis functions for adjoint source time function */
  src_fd = alloc1complexf(emf->nfreq);
  ntrace = 2*emf->nfreq;
  xx = alloc1float(emf->nt*ntrace);
  yy = alloc1float(2*emf->nfreq*ntrace);
  memset(xx, 0, emf->nt*ntrace*sizeof(float));
  for(j=0; j<ntrace; j++){
    for(ifreq=0; ifreq<2*emf->nfreq; ++ifreq) {
      i = ifreq + 2*emf->nfreq*j;
      yy[i] = (ifreq==j)?1.:0.;
    }
  }
  adjsrc_freq2time_init(emf->omega0, emf->dt, emf->omegas, emf->nfreq, emf->nt, ntrace, 
			adjsrc_tol, adjsrc_eps, adjsrc_nrestart, adjsrc_niter, emf->verb);
  adjsrc_init(xx, yy);
  adjsrc_freq2time();/* most computationally intensive part hidden here */    
  adjsrc_freq2time_close();

  /* 1. compute data misfit + adjoint source by linear combination of basis functions */
  fcost = 0.;
  memset(&emf->obj[0][0][0], 0, emf->nfreq*acqui->nrec*emf->nchrec*sizeof(float));
  memset(&emf->dres_td[0][0][0], 0, emf->nt*acqui->nrec*emf->nchrec*sizeof(float));
  for(ichrec=0; ichrec<emf->nchrec; ichrec++){
    for(irec=0; irec<acqui->nrec; irec++){
      if(emf->offset_ok[irec]){
	for(ifreq=0; ifreq<emf->nfreq; ++ifreq){
	  emf->dres_fd[ichrec][ifreq][irec] = emf->dobs_fd[ichrec][ifreq][irec] - emf->dcal_fd[ichrec][ifreq][irec];
	  emf->dres_fd[ichrec][ifreq][irec] /= emf->delta_emf[ichrec][ifreq][irec];/* W(dobs-R u) */
	  
	  emf->obj[ichrec][ifreq][irec] = cabs(emf->dres_fd[ichrec][ifreq][irec]);
	  s = emf->dres_fd[ichrec][ifreq][irec]*conj(emf->dres_fd[ichrec][ifreq][irec]);
	  src_fd[ifreq] = emf->dres_fd[ichrec][ifreq][irec]/emf->delta_emf[ichrec][ifreq][irec]; /* W^t*W(dobs-R u) */

	  /* do not use data smaller than noise floor */
	  if((strcmp(emf->chrec[ichrec],"Ex")==0||strcmp(emf->chrec[ichrec],"Ey")==0) &&
	     cabs(emf->dobs_fd[ichrec][ifreq][irec])<emf->noisefloorE){
	    src_fd[ifreq] = 0.;
	    s = 0.;
	    emf->obj[ichrec][ifreq][irec] = 0.;
	  }
	  if((strcmp(emf->chrec[ichrec],"Hx")==0||strcmp(emf->chrec[ichrec],"Hy")==0) &&
	     cabs(emf->dobs_fd[ichrec][ifreq][irec])<emf->noisefloorH){
	    src_fd[ifreq] = 0.;
	    s = 0.;
	    emf->obj[ichrec][ifreq][irec] = 0.;
	  }
	  fcost += s;
	}
	for(it=0; it<emf->nt; it++){
	  s = 0.;
	  for(ifreq=0; ifreq<emf->nfreq; ifreq++){
	    //minus sign due to conjugate of adjoint source
	    s += xx[it+emf->nt*ifreq]*crealf(src_fd[ifreq]);
	    s -= xx[it+emf->nt*(ifreq+emf->nfreq)]*cimagf(src_fd[ifreq]);
	  }
	  emf->dres_td[it][ichrec][irec] = s;//scale data residuals around 1
	}
	
      }
    }/* end for irec */
  }/* end for ichrec */
  fcost *= 0.5;

  free1float(xx);
  free1float(yy);
  free1complexf(src_fd);
  
  /* print out significant misfit */
  char fname[sizeof("sig_misfit_0000.txt")];
  sprintf(fname, "sig_misfit_%04d.txt", acqui->shot_idx[iproc]);
  write_misfit(acqui, emf, fname, emf->obj);
  
  return fcost; 
}

/*----------------------------------------------*/
void gn_adjoint_source(acqui_t *acqui, emf_t *emf)
/*< compute cost function and the adjoint source for FWI >*/
{
  int ichrec, irec, ifreq, it, ntrace, i, j;
  float _Complex s,omegap;
  float *xx, *yy;

  /* 0. compute basis functions for adjoint source time function */
  ntrace = 2*emf->nfreq;
  xx = alloc1float(emf->nt*ntrace);
  yy = alloc1float(2*emf->nfreq*ntrace);
  memset(xx, 0, emf->nt*ntrace*sizeof(float));
  for(j=0; j<ntrace; j++){
    for(ifreq=0; ifreq<2*emf->nfreq; ++ifreq) {
      i = ifreq + 2*emf->nfreq*j;
      yy[i] = (ifreq==j)?1.:0.;
    }
  }
  adjsrc_freq2time_init(emf->omega0, emf->dt, emf->omegas, emf->nfreq, emf->nt, ntrace, 
			adjsrc_tol, adjsrc_eps, adjsrc_nrestart, adjsrc_niter, emf->verb);
  adjsrc_init(xx, yy);
  adjsrc_freq2time();/* most computationally intensive part hidden here */    
  adjsrc_freq2time_close();

  //convert time domain w into frequency domain
  for(ichrec=0; ichrec<emf->nchrec; ++ichrec){
    for(irec = 0; irec<acqui->nrec; irec++) {
      for(ifreq=0; ifreq<emf->nfreq; ++ifreq){
	omegap = (1.0+I)*sqrt(emf->omega0*emf->omegas[ifreq]);
	s = 0;
	for(it=0; it<emf->nt; it++) s += emf->dres_td[it][ichrec][irec]*cexp(I*omegap*it*emf->dt);
	emf->dres_fd[ichrec][ifreq][irec] = s;
      }
    }
  }

  //apply frequency domain weighting
  memset(&emf->dres_td[0][0][0], 0, emf->nt*acqui->nrec*emf->nchrec*sizeof(float));
  for(ichrec=0; ichrec<emf->nchrec; ichrec++){
    for(irec=0; irec<acqui->nrec; irec++){
      if(emf->offset_ok[irec]){
	for(ifreq=0; ifreq<emf->nfreq; ++ifreq){
	  emf->dres_fd[ichrec][ifreq][irec] /= emf->delta_emf[ichrec][ifreq][irec];/* W(dobs-R u) */
	  s = emf->dres_fd[ichrec][ifreq][irec]/emf->delta_emf[ichrec][ifreq][irec]; /* W^t*W(dobs-R u) */

	  for(it=0; it<emf->nt; it++)
	    emf->dres_td[it][ichrec][irec] += xx[it+emf->nt*ifreq]*creal(s) - xx[it+emf->nt*(ifreq+emf->nfreq)]*cimag(s);
	}
      }
    }/* end for irec */
  }/* end for ichrec */

  free1float(xx);
  free1float(yy);

}
