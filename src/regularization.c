/* regularization by Tikhonov or total variation (TV)
 *------------------------------------------------------------------------
 *
 * Copyright (c) 2021 Harbin Institute of Technology. All rights reserved.
 * Anothr: Pengliang Yang 
 * Email: ypl.2100@gmail.com
 * Homepage: https://yangpl.wordpress.com
 *-----------------------------------------------------------------------*/
#include "cstd.h"

/*< c_h and c_v for regularization in horizontal and vertical directions >*/
float regularization_tikhonov(float *x, float *g, int n1, int n2, int n3, float d1, float d2, float d3)
{
  int i1, i2, i3, k0, kp1, km1;
  float s1, s2, s3, t1, t2, t3, tmp, fcost_mod;
  float _d1 = 1.;///(d1*d1);
  float _d2 = 1.;///(d2*d2);
  float _d3 = 1.;///(d3*d3);
  float c_h = 1;//large coefficient to penalize horizontal changes
  float c_v = 0.03;//small coefficients to panalize vertical changes

  memset(g, 0, n1*n2*n3*sizeof(float));
  
  s1 = 0.;
  s2 = 0.;
  s3 = 0.;
  for(i3=1; i3<n3-1; i3++){
    for(i2=1; i2<n2-1; i2++){
      for(i1=1; i1<n1-1; i1++){
	k0 = i1 + n1*(i2 + n2*i3);

	km1 = k0-1;
	kp1 = k0+1;
	t1 = -(x[km1] -2.0*x[k0]+x[kp1])*_d1;
	tmp = x[kp1] - x[k0];
	s1 += tmp*tmp;

	km1 = k0-n1;
	kp1 = k0+n1;
	t2 = -(x[km1] -2.0*x[k0]+x[kp1])*_d2;
	tmp = x[kp1] - x[k0];
	s2 += tmp*tmp;

	km1 = k0-n1*n2;
	kp1 = k0+n1*n2;
	t3 = -(x[km1] -2.0*x[k0]+x[kp1])*_d3;
	tmp = x[kp1] - x[k0];
	s3 += tmp*tmp;

	g[k0] = c_h*(t1+t2)+c_v*t3;
      }
    }
  }
  s1 *= _d1;
  s2 *= _d2;
  s3 *= _d3;
  
  fcost_mod = c_h*(s1 + s2) + c_v*s3;

  return fcost_mod;
}


float regularization_tv(float *x, float *g, int n1, int n2, int n3, float d1, float d2, float d3)
{
  int i1, i2, i3, k0, kp1, km1;
  float s1, s2, s3, tmp, fcost_mod;
  float _d1 = 1.; //(d1*d1);
  float _d2 = 1.; //(d2*d2);
  float _d3 = 1.;///(d3*d3);
  float beta = 1e-14;

  float *g1 = alloc1float(n1*n2*n3);
  float *g2 = alloc1float(n1*n2*n3);
  float *g3 = alloc1float(n1*n2*n3);

  memset(g, 0, n1*n2*n3*sizeof(float));
  memset(g1, 0, n1*n2*n3*sizeof(float));
  memset(g2, 0, n1*n2*n3*sizeof(float));
  memset(g3, 0, n1*n2*n3*sizeof(float));

  fcost_mod = 0;
  for(i3=0; i3<n3-1; i3++){
    for(i2=0; i2<n2-1; i2++){
      for(i1=0; i1<n1-1; i1++){
	k0 = i1 + n1*(i2 + n2*i3);

	kp1 = k0 + 1;
	s1 = (x[kp1] - x[k0])*_d1;//du/dx

	kp1 = k0 + n1;
	s2 = (x[kp1] - x[k0])*_d2;//du/dy

	kp1 = k0 + n1*n2;
	s3 = (x[kp1] - x[k0])*_d3;//du/dz

	tmp = sqrt(s1*s1 + s2*s2 + s3*s3 + beta);//sqrt(|grad|^2 + beta)
	fcost_mod += tmp;
	//grad/sqrt(|grad|^2 + beta)
	g1[k0] = s1/tmp;
	g2[k0] = s2/tmp;
	g3[k0] = s3/tmp;
      }
    }
  }
  
  for(i3=1; i3<n3; i3++){
    for(i2=1; i2<n2; i2++){
      for(i1=1; i1<n1; i1++){
	k0 = i1 + n1*(i2 + n2*i3);

	km1 = k0 - 1;
	s1 = (g1[k0] - g1[km1])*_d1;

	km1 = k0 - n1;
	s2 = (g2[k0] - g2[km1])*_d2;

	km1 = k0 - n1*n2;
	s3 = (g3[k0] - g3[km1])*_d3;

	g[k0] = -(s1 + s2 + s3);  //-div(grad/|grad|)
      }
    }
  }

  free1float(g1);
  free1float(g2);
  free1float(g3);
  
  return fcost_mod;
}
  
