/* build inversion gradient for EM inversion
 *------------------------------------------------------------------------
 *
 * Copyright (c) 2021 Harbin Institute of Technology. All rights reserved.
 * Anothr: Pengliang Yang 
 * Email: ypl.2100@gmail.com
 * Homepage: https://yangpl.wordpress.com
 *-----------------------------------------------------------------------*/
#include "cstd.h"
#include "emf.h"
#include "fwi.h"

void build_gradient(emf_t *emf, fwi_t *fwi)
{
  int ifreq, i1, i2, i3, j1, j2, j3;
  float _Complex fwd_emf, adj_emf;
  float gh, gv, d3s, d3n;
		    
  for(i3=0; i3<emf->n3; ++i3){
    j3 = i3+ emf->nbe;
    d3n = emf->nugrid?(emf->x3s[j3]-emf->x3s[j3-1]):emf->d3;//centered at nonstaggered grid
    d3s = emf->nugrid?(emf->x3n[j3+1]-emf->x3n[j3]):emf->d3;//centered at staggered grid
    for(i2=0; i2<emf->n2; ++i2){
      j2 = i2+ emf->nbe;
      for(i1=0; i1<emf->n1; ++i1){
  	j1 = i1+ emf->nbe;

  	gh = 0.0;
	gv = 0.0;
  	for(ifreq=0; ifreq<emf->nfreq; ++ifreq){
	  /*-------------------------------------------------------*/
	  if(emf->rd1==1){
	    fwd_emf = 0.5*(emf->fwd_E1[ifreq][j3][j2][j1] +emf->fwd_E1[ifreq][j3][j2][j1-1]);
	    adj_emf = 0.5*(emf->adj_E1[ifreq][j3][j2][j1] +emf->adj_E1[ifreq][j3][j2][j1-1]);
	  }else if(emf->rd1==2){
	    /* we use 4-th order averaging on staggered grid: c1=9/16, c2=-1/16 */
	    /* f(x)=c1*(f(x+0.5*h)-f(x-0.5*h)) + c2*(f(x+1.5*h)-f(x-1.5*h))  */
	    fwd_emf = 0.5625*(emf->fwd_E1[ifreq][j3][j2][j1] +emf->fwd_E1[ifreq][j3][j2][j1-1])
	      - 0.0625*(emf->fwd_E1[ifreq][j3][j2][j1+1] +emf->fwd_E1[ifreq][j3][j2][j1-2]);
	    adj_emf = 0.5625*(emf->adj_E1[ifreq][j3][j2][j1] +emf->adj_E1[ifreq][j3][j2][j1-1])
	      - 0.0625*(emf->adj_E1[ifreq][j3][j2][j1+1] +emf->adj_E1[ifreq][j3][j2][j1-2]);
	  }
  	  gh += creal(fwd_emf*adj_emf);
	  /*-------------------------------------------------------*/
	  if(emf->rd2==1){
	    fwd_emf = 0.5*(emf->fwd_E2[ifreq][j3][j2][j1] +emf->fwd_E2[ifreq][j3][j2-1][j1]);
	    adj_emf = 0.5*(emf->adj_E2[ifreq][j3][j2][j1] +emf->adj_E2[ifreq][j3][j2-1][j1]);
	  }else if(emf->rd2==2){
	    fwd_emf = 0.5625*(emf->fwd_E2[ifreq][j3][j2][j1] +emf->fwd_E2[ifreq][j3][j2-1][j1])
	      - 0.0625*(emf->fwd_E2[ifreq][j3][j2+1][j1] +emf->fwd_E2[ifreq][j3][j2-2][j1]);
	    adj_emf = 0.5625*(emf->adj_E2[ifreq][j3][j2][j1] +emf->adj_E2[ifreq][j3][j2-1][j1])
	      - 0.0625*(emf->adj_E2[ifreq][j3][j2+1][j1] +emf->adj_E2[ifreq][j3][j2-2][j1]);
	  }
  	  gh += creal(fwd_emf*adj_emf);
	  /*-------------------------------------------------------*/
	  if(emf->rd3==1){
	    fwd_emf = 0.5*(  emf->fwd_E3[ifreq][j3][j2][j1]/emf->inveps33[j3][j2][j1]
	  		     + emf->fwd_E3[ifreq][j3-1][j2][j1]/emf->inveps33[j3-1][j2][j1]);
	    adj_emf = 0.5*(  emf->adj_E3[ifreq][j3][j2][j1]/emf->inveps33[j3][j2][j1]
	  		     + emf->adj_E3[ifreq][j3-1][j2][j1]/emf->inveps33[j3-1][j2][j1]);
	  }else if(emf->rd3==2){
	    /* Ez is discontinous in z, but Jz is continuous, average over Jz and then back to Ez */
	    fwd_emf = 0.5625*(  emf->fwd_E3[ifreq][j3][j2][j1]/emf->inveps33[j3][j2][j1]
	  			+ emf->fwd_E3[ifreq][j3-1][j2][j1]/emf->inveps33[j3-1][j2][j1])
	      - 0.0625*(  emf->fwd_E3[ifreq][j3+1][j2][j1]/emf->inveps33[j3+1][j2][j1]
	  		  + emf->fwd_E3[ifreq][j3-2][j2][j1]/emf->inveps33[j3-2][j2][j1]);
	    adj_emf = 0.5625*(  emf->adj_E3[ifreq][j3][j2][j1]/emf->inveps33[j3][j2][j1]
	  			+ emf->adj_E3[ifreq][j3-1][j2][j1]/emf->inveps33[j3-1][j2][j1])
	      - 0.0625*(  emf->adj_E3[ifreq][j3+1][j2][j1]/emf->inveps33[j3+1][j2][j1]
	  		  + emf->adj_E3[ifreq][j3-2][j2][j1]/emf->inveps33[j3-2][j2][j1]);
	  }
	  fwd_emf *= 0.5*(emf->inveps33[j3][j2][j1] + emf->inveps33[j3-1][j2][j1]);
	  adj_emf *= 0.5*(emf->inveps33[j3][j2][j1] + emf->inveps33[j3-1][j2][j1]);
	  gv += creal(fwd_emf*adj_emf);
	  /*-------------------------------------------------------*/
  	}
	
  	/* multiply the volume such that the gradient has correct magnitude */
	fwi->g_h[i3][i2][i1] = gh*emf->d1*emf->d2*d3n/emf->rho_h[i3][i2][i1];//dJ/d(log(rho))=<E^a|E^f>/rho
	fwi->g_v[i3][i2][i1] = gv*emf->d1*emf->d2*d3n/emf->rho_v[i3][i2][i1];//dJ/d(log(rho))=<E^a|E^f>/rho

	if(i3<emf->ibathy[i2][i1]){
	  fwi->g_h[i3][i2][i1] = 0.;
	  fwi->g_v[i3][i2][i1] = 0.;
	}
      }
    }
  }

}
