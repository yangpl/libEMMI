/* Compute effective medium by homogenization/averaging by Davedycheva 2003
 *
 * Key point of the method:
 *     resistivity connected in series vertically, connected in parallel horizontally
 * This leads to the following averaging scheme:
 *    1) average conductivity in horizontal direction 
 *    2) average resistivity in vertical direction, 
 *
 *------------------------------------------------------------------------------
 *
 * Copyright (c) 2021 Harbin Institute of Technology. All rights reserved.
 * Anothr: Pengliang Yang 
 * Email: ypl.2100@gmail.com
 * Homepage: https://yangpl.wordpress.com
 *-----------------------------------------------------------------------------*/
#include "cstd.h"
#include "emf.h"

void homogenization(emf_t *emf)
{
  int i1, i2, i3;
  float tmp;

  //do homogenization for the resistivity below bathymetry
  for(i3=0; i3<emf->n3; i3++){
    for(i2=0; i2<emf->n2; i2++){
      for(i1=0; i1<emf->n1-1; i1++){
	if(i3>emf->ibathy[i2][i1]){
	  tmp = 0.5*(1./emf->rho[i3][i2][i1] + 1./emf->rho[i3][i2][i1+1]);
	  emf->rho11[i3][i2][i1] = 1./tmp;
	}
      }
      i1 = emf->n1-1;
      if(i3>emf->ibathy[i2][i1]) emf->rho11[i3][i2][i1] = emf->rho[i3][i2][i1];
    }
  }

  for(i3=0; i3<emf->n3; i3++){
    for(i1=0; i1<emf->n1; i1++){
      for(i2=0; i2<emf->n2-1; i2++){
	if(i3>emf->ibathy[i2][i1]){
	  tmp = 0.5*(1./emf->rho[i3][i2][i1] + 1./emf->rho[i3][i2+1][i1]);
	  emf->rho22[i3][i2][i1] = 1./tmp;
	}
      }
      i2 = emf->n2-1;
      if(i3>emf->ibathy[i2][i1]) emf->rho22[i3][i2][i1] = emf->rho[i3][i2][i1];
    }
  }
  
  for(i2=0; i2<emf->n2; i2++){
    for(i1=0; i1<emf->n1; i1++){
      for(i3=emf->ibathy[i2][i1]+1; i3<emf->n3-1; i3++){
	if(emf->nugrid){
	  float s0 = emf->x3n[i3+1] - emf->x3n[i3];
	  float s1 = emf->x3s[i3] - emf->x3n[i3];
	  float s2 = emf->x3n[i3+1] - emf->x3s[i3];
	  emf->rho33[i3][i2][i1] = (s1*emf->rho[i3][i2][i1] + s2*emf->rho[i3+1][i2][i1])/s0;
	}else
	  emf->rho33[i3][i2][i1] = 0.5*(emf->rho[i3][i2][i1] + emf->rho[i3+1][i2][i1]);
      }
      emf->rho33[emf->n3-1][i2][i1] = emf->rho[emf->n3-1][i2][i1];
    }
  }
  
}

