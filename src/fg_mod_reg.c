/* add model regularization term for function and gradient evaluation
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
#include "constants.h"

float regularization_tikhonov(float *x, float *g, int n1, int n2, int n3, float d1, float d2, float d3);
float regularization_tv(float *x, float *g, int n1, int n2, int n3, float d1, float d2, float d3);

/* Harlan 1995, regularization by model reparameterization, 
 * We should always use gradient smoothing to achieve faster convergence **/
void triangle_smoothing(float ***mod, int n1, int n2, int n3, int r1, int r2, int r3, int repeat);


/*---------------------------------------------------------------*/
void fg_mod_reg(emf_t *emf, fwi_t *fwi, float *x, float *g)
/*< model regularization (assume fwi gradient has been computed and stored in grad) >*/
{
  FILE *fp;
  int i1, i2, i3, i, j, k;
  float gamma1, gamma2;
    
  float beta = pow(0.9, fwi->iter);//beta is a cooling factor
  float fcost_mod1 = 0;
  float fcost_mod2 = 0;
  float *g1 = alloc1float(fwi->n);
  float *g2 = alloc1float(fwi->n);
  float *xdiff = alloc1float(fwi->n);
  beta *= fwi->ndp;
  
  /* /\* step 1: gradient preconditioning by smoothing *\/ */
  /* /\* mirror around bathymetry before smoothing, keep the average around bathymetry *\/ */
  /* memcpy(&fwi->grad[0][0][0], g, emf->n123*sizeof(float)); */
  /* for(i2 = 0; i2<emf->n2; i2++){ */
  /*   for(i1 = 0; i1<emf->n1; i1++){ */
  /*     for(j=0,i3=emf->ibathy[i2][i1]-1; i3>=0; i3--,j++){ */
  /* 	k = MIN(emf->ibathy[i2][i1] + j, emf->n3-1); */
  /* 	fwi->grad[i3][i2][i1] = fwi->grad[k][i2][i1]; */
  /*     } */
  /*   } */
  /* } */
  /* triangle_smoothing(fwi->grad, emf->n1, emf->n2, emf->n3, fwi->r1, fwi->r2, fwi->r3, fwi->repeat); */
  /* memcpy(g, &fwi->grad[0][0][0], emf->n123*sizeof(float)); */

  
  /* step 2: Tikhonov regularization + TV regularization */
  fwi->fcost_mod = 0;
  memset(g1, 0, fwi->n*sizeof(float));
  memset(g2, 0, fwi->n*sizeof(float));
  for(i3=0; i3<emf->n3; i3++){
    for(i2=0; i2<emf->n2; i2++){
      for(i1=0; i1<emf->n1; i1++){
	i = i1 + emf->n1*(i2 + emf->n2*i3);
	xdiff[i] = x[i] - fwi->xref[i];
      }
    }
  }
  fcost_mod1 = regularization_tikhonov(xdiff, g1, emf->n1, emf->n2, emf->n3, emf->d1, emf->d2, emf->d3);
  fcost_mod2 = regularization_tv(xdiff, g2, emf->n1, emf->n2, emf->n3, emf->d1, emf->d2, emf->d3);

  /* step 3: add regularization temse in fcost and gradient */
  if(fwi->gamma1>0) {
    gamma1 = fwi->gamma1*beta;
    fcost_mod1 *= gamma1;
    fwi->fcost_mod += fcost_mod1;
    for(i2 = 0; i2<emf->n2; i2++){
      for(i1 = 0; i1<emf->n1; i1++){
	for(i3=emf->ibathy[i2][i1]+1; i3<emf->n3; i3++){
	  i = i1 + emf->n1*(i2 + emf->n3*i3);
	  g[i] += g1[i]*gamma1;
	}
      }
    }
    
    if(emf->verb){
      printf("fcost_Tikhonov=%g\n", fcost_mod1);

      fp = fopen("gradient_tikhonov", "wb");
      fwrite(g1, emf->n123*sizeof(float), 1, fp);
      fclose(fp);
    }
  }
  if(fwi->gamma2>0) {
    gamma2 = fwi->gamma2*beta;
    fcost_mod2 *= gamma2;
    fwi->fcost_mod += fcost_mod2;
    for(i2 = 0; i2<emf->n2; i2++){
      for(i1 = 0; i1<emf->n1; i1++){
	for(i3=emf->ibathy[i2][i1]+1; i3<emf->n3; i3++){
	  i = i1 + emf->n1*(i2 + emf->n3*i3);
	  g[i] += g2[i]*gamma2;
	}
      }
    }
    if(emf->verb){
      printf("fcost_TV=%g\n", fcost_mod2);

      fp = fopen("gradient_tv", "wb");
      fwrite(g2, emf->n123*sizeof(float), 1, fp);
      fclose(fp);
    }
  }

  for(i3=0; i3<emf->n3; i3++){
    for(i2=0; i2<emf->n2; i2++){
      for(i1=0; i1<emf->n1; i1++){
	k = i1 + emf->n1*(i2 + emf->n2*i3);
	//do not update water column and the 1st formation layer
	if(i3<= emf->ibathy[i2][i1]) g[k] = 0.;
	if(emf->invmask && !(emf->mask[i3][i2][i1]>0)) g[k] = 0;
      }
    }
  }

  free1float(xdiff);
  free1float(g1);
  free1float(g2);

}


/*---------------------------------------------------------------*/
void Hv_mod_reg(emf_t *emf, fwi_t *fwi, float *r, float *Hv)
/*< model regularization (assume fwi gradient has been computed and stored in grad) >*/
{
  int i1, i2, i3, k0, kp1, km1;
  float gamma1, t1, t2, t3;
    
  float beta = pow(0.9, fwi->iter);//beta is a cooling factor
  float _d1 = 1.;///(d1*d1);
  float _d2 = 1.;///(d2*d2);
  float _d3 = 1.;///(d3*d3);
  float c_h = 1;//large coefficient to penalize horizontal changes
  float c_v = 0.03;//small coefficients to panalize vertical changes
  //beta *= fwi->ndp;
  gamma1 = fwi->gamma1*beta;
  
  for(i3=1; i3<emf->n3-1; i3++){
    for(i2=1; i2<emf->n2-1; i2++){
      for(i1=1; i1<emf->n1-1; i1++){

	k0 = i1 + emf->n1*(i2 + emf->n2*i3);
	if(i3>emf->ibathy[i2][i1]){
	  km1 = k0-1;
	  kp1 = k0+1;
	  t1 = -(r[km1] -2.0*r[k0]+r[kp1])*_d1;

	  km1 = k0-emf->n1;
	  kp1 = k0+emf->n1;
	  t2 = -(r[km1] -2.0*r[k0]+r[kp1])*_d2;

	  km1 = k0-emf->n1*emf->n2;
	  kp1 = k0+emf->n1*emf->n2;
	  t3 = -(r[km1] -2.0*r[k0]+r[kp1])*_d3;

	  Hv[k0] += (c_h*(t1+t2)+c_v*t3)*gamma1;
	}
      }
    }
  }

}
