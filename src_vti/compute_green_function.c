/* compute Green's function Gij=Ei/Jj
 *------------------------------------------------------------------------
 *
 * Copyright (c) 2021 Harbin Institute of Technology. All rights reserved.
 * Anothr: Pengliang Yang 
 * Email: ypl.2100@gmail.com
 *-----------------------------------------------------------------------*/
#include "cstd.h"
#include "emf.h"

void compute_green_function(emf_t *emf, int adj)
{
  int ifreq,it,i1,i2,i3;
  float _Complex s, omegap, src_fd;

  int i1min=emf->nb;
  int i2min=emf->nb;
  int i3min=emf->nb;
  int i1max=emf->n1pad-1-emf->nb;
  int i2max=emf->n2pad-1-emf->nb;
  int i3max=emf->n3pad-1-emf->nb;

  for(ifreq=0; ifreq < emf->nfreq; ++ifreq){
    /*------------------- DTFT over omega'-------------------------*/
    omegap = (1.0+I)*sqrtf(emf->omega0*emf->omegas[ifreq]);/* omega' in fictitous wave domain */
    src_fd = 0.;
    for(it=0; it<emf->nt; ++it) src_fd += emf->stf[it]*cexp(I*omegap*it*emf->dt);//J'
    s = csqrt(-I*0.5*emf->omegas[ifreq]/emf->omega0);

    if(adj){
      for(i3=i3min; i3<=i3max; ++i3){
	for(i2=i2min; i2<=i2max; ++i2){
	  for(i1=i1min; i1<=i1max; ++i1){
	    emf->adj_E1[ifreq][i3][i2][i1] /= src_fd;
	    emf->adj_E2[ifreq][i3][i2][i1] /= src_fd;
	    emf->adj_E3[ifreq][i3][i2][i1] /= src_fd;
	    emf->adj_E1[ifreq][i3][i2][i1] *= s;
	    emf->adj_E2[ifreq][i3][i2][i1] *= s;
	    emf->adj_E3[ifreq][i3][i2][i1] *= s;
	  }
	}
      }

    }else{
      for(i3=i3min; i3<=i3max; ++i3){
	for(i2=i2min; i2<=i2max; ++i2){
	  for(i1=i1min; i1<=i1max; ++i1){
	    emf->fwd_E1[ifreq][i3][i2][i1] /= src_fd;
	    emf->fwd_E2[ifreq][i3][i2][i1] /= src_fd;
	    emf->fwd_E3[ifreq][i3][i2][i1] /= src_fd;
	    emf->fwd_E1[ifreq][i3][i2][i1] *= s;
	    emf->fwd_E2[ifreq][i3][i2][i1] *= s;
	    emf->fwd_E3[ifreq][i3][i2][i1] *= s;
	  }
	}
      }

    }

  }
}

