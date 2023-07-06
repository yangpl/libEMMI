/* Discrete time Fourier transform (DTFT) of EM fields using wave-diffusion transform
 *------------------------------------------------------------------------
 *
 * Copyright (c) 2021 Harbin Institute of Technology. All rights reserved.
 * Anothr: Pengliang Yang 
 * Email: ypl.2100@gmail.com
 * Homepage: https://yangpl.wordpress.com
 *-----------------------------------------------------------------------*/
#include "cstd.h"
#include "emf.h"
#include "constants.h"
#include "mpi_info.h"

/*-------------------------------------------------------------*/
void dtft_emf_init(emf_t *emf, int adj)
{
  if(adj){
    emf->adj_E1 = alloc4complexf(emf->n1pad, emf->n2pad, emf->n3pad, emf->nfreq);
    emf->adj_E2 = alloc4complexf(emf->n1pad, emf->n2pad, emf->n3pad, emf->nfreq);
    emf->adj_E3 = alloc4complexf(emf->n1pad, emf->n2pad, emf->n3pad, emf->nfreq);
    memset(&emf->adj_E1[0][0][0][0], 0, emf->nfreq*emf->n123pad*sizeof(float _Complex));
    memset(&emf->adj_E2[0][0][0][0], 0, emf->nfreq*emf->n123pad*sizeof(float _Complex));
    memset(&emf->adj_E3[0][0][0][0], 0, emf->nfreq*emf->n123pad*sizeof(float _Complex));
  }else{
    emf->fwd_E1 = alloc4complexf(emf->n1pad, emf->n2pad, emf->n3pad, emf->nfreq);
    emf->fwd_E2 = alloc4complexf(emf->n1pad, emf->n2pad, emf->n3pad, emf->nfreq);
    emf->fwd_E3 = alloc4complexf(emf->n1pad, emf->n2pad, emf->n3pad, emf->nfreq);
    memset(&emf->fwd_E1[0][0][0][0], 0, emf->nfreq*emf->n123pad*sizeof(float _Complex));
    memset(&emf->fwd_E2[0][0][0][0], 0, emf->nfreq*emf->n123pad*sizeof(float _Complex));
    memset(&emf->fwd_E3[0][0][0][0], 0, emf->nfreq*emf->n123pad*sizeof(float _Complex));
  }
}

void dtft_emf_close(emf_t *emf, int adj)
{
  if(adj){
    free4complexf(emf->adj_E1);
    free4complexf(emf->adj_E2);
    free4complexf(emf->adj_E3);
  }else{
    free4complexf(emf->fwd_E1);
    free4complexf(emf->fwd_E2);
    free4complexf(emf->fwd_E3);
  }

}

/*-----------------------------------------------------------*/
void dtft_emf(emf_t *emf, int it, int adj)
{
  int i1, i2, i3, ifreq;
  float _Complex factor, omegap;
  
  int i1min = emf->nb;
  int i2min = emf->nb;
  int i3min = emf->nb;
  int i1max = emf->n1pad-1-emf->nb;
  int i2max = emf->n2pad-1-emf->nb;
  int i3max = emf->n3pad-1-emf->nb;

  for(ifreq=0; ifreq<emf->nfreq; ++ifreq){
    omegap = (1.0+I)*sqrt(emf->omega0*emf->omegas[ifreq]);
    factor = cexp(I*omegap*(it+0.5)*emf->dt);
    
    if(adj){
#ifdef _OPENMP
#pragma omp parallel for default(none)					\
  schedule(static)							\
  private(i1, i2, i3)							\
  shared(i1min, i1max, i2min, i2max, i3min, i3max, ifreq, factor,emf)
#endif
      for(i3=i3min; i3<=i3max; ++i3){
	for(i2=i2min; i2<=i2max; ++i2){
	  for(i1=i1min; i1<=i1max; ++i1){
	    emf->adj_E1[ifreq][i3][i2][i1] += emf->E1[i3][i2][i1]*factor;
	    emf->adj_E2[ifreq][i3][i2][i1] += emf->E2[i3][i2][i1]*factor;
	    emf->adj_E3[ifreq][i3][i2][i1] += emf->E3[i3][i2][i1]*factor;
	  }
	}
      }


    }else{
#ifdef _OPENMP
#pragma omp parallel for default(none)					\
  schedule(static)							\
  private(i1, i2, i3)							\
  shared(i1min, i1max, i2min, i2max, i3min, i3max, ifreq, factor,emf)
#endif
      for(i3=i3min; i3<=i3max; ++i3){
	for(i2=i2min; i2<=i2max; ++i2){
	  for(i1=i1min; i1<=i1max; ++i1){
	    emf->fwd_E1[ifreq][i3][i2][i1] += emf->E1[i3][i2][i1]*factor;
	    emf->fwd_E2[ifreq][i3][i2][i1] += emf->E2[i3][i2][i1]*factor;
	    emf->fwd_E3[ifreq][i3][i2][i1] += emf->E3[i3][i2][i1]*factor;
	  }
	}
      }

    }
  }/* end for ifreq */
}
