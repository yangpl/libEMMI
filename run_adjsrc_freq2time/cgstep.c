/* Claerbout's conjugate-gradient iteration */
/*
  Copyright (C) 2004 University of Texas at Austin
  
  This program is free software; you can redistribute it and/or modify
  it under the temse of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.
  
  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
  
  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/
#include "cstd.h"

static const double EPSILON=1.e-12; /* precision */

static float * S;  /* model step */
static float * Ss; /* residual step */
static bool Allocated = false; /* if S and Ss are allocated */

static float  dotprod (int n, const float * x, const float * y)
/* float  dot product */
{
  float  prod;
  int i;

  prod = 0.0;
  for (i = 0; i < n; i++) prod += x[i]*y[i];
  return prod;
}

void sf_cgstep( bool forget             /* restart flag */, 
		 int nx                  /* model size */, 
		 int ny                  /* data size */, 
		 float * x        /* current model [nx] */,  
		 const float * g  /* gradient [nx] */, 
		 float * rr       /* data residual [ny] */,
		 const float * gg /* conjugate gradient [ny] */) 
/*< Step of Claerbout's conjugate-gradient iteration for float  operators. 
  The data residual is rr = A x - dat  >*/
{
  float  sds, gdg, gds, sdg, determ, gdr, sdr, alfa, beta;
  int i;
  if (Allocated == false) {
    Allocated = forget = true;
    S  = alloc1float (nx);
    Ss = alloc1float (ny);
  }
  if (forget) {
    for (i = 0; i < nx; i++) {
      S[i] = 0.0;
    }
    for (i = 0; i < ny; i++) {
      Ss[i] = 0.0;
    }    

    beta = 0.0;
    if (dotprod(ny, gg, gg) <= 0.) return;
    alfa = - dotprod( ny, gg, rr) / dotprod(ny, gg, gg);
  } else {
    /* search plane by solving 2-by-2
       G . (R - G*alfa - S*beta) = 0
       S . (R - G*alfa - S*beta) = 0 */
    gdg = dotprod( ny, gg, gg);       
    sds = dotprod( ny, Ss, Ss);       
    gds = dotprod( ny, gg, Ss);    
    sdg = dotprod( ny, Ss, gg);   
    if (fabs(gdg) == 0. || fabs(sds) == 0.) return;

    determ = 1.0 - (gds/gdg)*(gds/sds);
    if (creal(determ) > EPSILON) {
      determ *= gdg * sds;
    } else {
      determ = gdg * sds * EPSILON;
    }

    gdr = - dotprod( ny, gg, rr);
    sdr = - dotprod( ny, Ss, rr);
    alfa = ( sds * gdr - gds * sdr ) / determ;
    beta = (-sdg * gdr + gdg * sdr ) / determ;
  }

  for (i = 0; i < nx; i++) {
    S[i]  =  alfa * g[i] + beta *  S[i];
    x[i] +=  S[i];
  }
  for (i = 0; i < ny; i++) {
    Ss[i] = alfa * gg[i] + beta * Ss[i];
    rr[i] += Ss[i];
  }
}

void sf_cgstep_close (void) 
/*< Free allocated space. >*/ 
{
  if (Allocated) {
    free (S);
    free (Ss);
    Allocated = false;
  }
}
