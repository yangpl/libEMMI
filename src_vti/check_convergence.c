/* check convergence of frequency domain EM fields at 8 corners of the cube
 *------------------------------------------------------------------------
 *
 * Copyright (c) 2021 Harbin Institute of Technology. All rights reserved.
 * Anothr: Pengliang Yang 
 * Email: ypl.2100@gmail.com
 *-----------------------------------------------------------------------*/
#include "cstd.h"
#include "emf.h"

int check_convergence(emf_t *emf, int adj)
{
  static float _Complex old1[8], new1[8], old2[8], new2[8];
  static int first = 1;
  int j, k;
  float _Complex ****E1, ****E2;
  
  if(adj){
    E1 = emf->adj_E1;
    E2 = emf->adj_E2;
  }else{
    E1 = emf->fwd_E1;
    E2 = emf->fwd_E2;
  }
  
  k = 0;
  if(first){
    old1[0] = E1[0][emf->nbe][emf->nbe][emf->nbe];
    old1[1] = E1[0][emf->nbe][emf->nbe][emf->nbe+emf->n1-1];
    old1[2] = E1[0][emf->nbe][emf->nbe+emf->n2-1][emf->nbe];
    old1[3] = E1[0][emf->nbe][emf->nbe+emf->n2-1][emf->nbe+emf->n1-1];
    old1[4] = E1[0][emf->nbe+emf->n3-1][emf->nbe][emf->nbe];
    old1[5] = E1[0][emf->nbe+emf->n3-1][emf->nbe][emf->nbe+emf->n1-1];
    old1[6] = E1[0][emf->nbe+emf->n3-1][emf->nbe+emf->n2-1][emf->nbe];
    old1[7] = E1[0][emf->nbe+emf->n3-1][emf->nbe+emf->n2-1][emf->nbe+emf->n1-1];
    
    old2[0] = E2[0][emf->nbe][emf->nbe][emf->nbe];
    old2[1] = E2[0][emf->nbe][emf->nbe][emf->nbe+emf->n1-1];
    old2[2] = E2[0][emf->nbe][emf->nbe+emf->n2-1][emf->nbe];
    old2[3] = E2[0][emf->nbe][emf->nbe+emf->n2-1][emf->nbe+emf->n1-1];
    old2[4] = E2[0][emf->nbe+emf->n3-1][emf->nbe][emf->nbe];
    old2[5] = E2[0][emf->nbe+emf->n3-1][emf->nbe][emf->nbe+emf->n1-1];
    old2[6] = E2[0][emf->nbe+emf->n3-1][emf->nbe+emf->n2-1][emf->nbe];
    old2[7] = E2[0][emf->nbe+emf->n3-1][emf->nbe+emf->n2-1][emf->nbe+emf->n1-1];
    
    first = 0;
  }else{
    new1[0] = E1[0][emf->nbe][emf->nbe][emf->nbe];
    new1[1] = E1[0][emf->nbe][emf->nbe][emf->nbe+emf->n1-1];
    new1[2] = E1[0][emf->nbe][emf->nbe+emf->n2-1][emf->nbe];
    new1[3] = E1[0][emf->nbe][emf->nbe+emf->n2-1][emf->nbe+emf->n1-1];
    new1[4] = E1[0][emf->nbe+emf->n3-1][emf->nbe][emf->nbe];
    new1[5] = E1[0][emf->nbe+emf->n3-1][emf->nbe][emf->nbe+emf->n1-1];
    new1[6] = E1[0][emf->nbe+emf->n3-1][emf->nbe+emf->n2-1][emf->nbe];
    new1[7] = E1[0][emf->nbe+emf->n3-1][emf->nbe+emf->n2-1][emf->nbe+emf->n1-1];
    
    new2[0] = E2[0][emf->nbe][emf->nbe][emf->nbe];
    new2[1] = E2[0][emf->nbe][emf->nbe][emf->nbe+emf->n1-1];
    new2[2] = E2[0][emf->nbe][emf->nbe+emf->n2-1][emf->nbe];
    new2[3] = E2[0][emf->nbe][emf->nbe+emf->n2-1][emf->nbe+emf->n1-1];
    new2[4] = E2[0][emf->nbe+emf->n3-1][emf->nbe][emf->nbe];
    new2[5] = E2[0][emf->nbe+emf->n3-1][emf->nbe][emf->nbe+emf->n1-1];
    new2[6] = E2[0][emf->nbe+emf->n3-1][emf->nbe+emf->n2-1][emf->nbe];
    new2[7] = E2[0][emf->nbe+emf->n3-1][emf->nbe+emf->n2-1][emf->nbe+emf->n1-1];

    for(j=0; j<8; j++){
      if(cabs(new1[j]-old1[j])<1e-2*cabs(new1[j]) && cabs(new2[j]-old2[j])<1e-2*cabs(new2[j])) k++;
      old1[j] = new1[j];
      old2[j] = new2[j];
    }
  }
  return k;
}

