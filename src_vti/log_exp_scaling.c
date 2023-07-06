/* log-exp switch between inversion parameter and physical parameter 
 *------------------------------------------------------------------------
 *
 * Copyright (c) 2021 Harbin Institute of Technology. All rights reserved.
 * Anothr: Pengliang Yang 
 * Email: ypl.2100@gmail.com
 * Homepage: https://yangpl.wordpress.com
 *-----------------------------------------------------------------------*/
#include <math.h>

void log_exp_scaling(int n, float *x, float *y, int flag)
{
  int i;
  if(flag==1){/*log transform from physical value to dimensionless value*/
    for(i=0; i<n; i++) y[i]=logf(x[i]);	
  }else if(flag==2){/*exp transform from dimensionless value to physical value*/
    for(i=0; i<n; i++) y[i]=expf(x[i]);	
  }
}
