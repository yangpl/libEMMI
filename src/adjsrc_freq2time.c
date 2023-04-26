/* convert adjoint source from freq to time domain via least-squares minimization: 
 *         min_m |Lm - d|^2 +eps^2*|Dm|  (implies:  Lm ~= d; Dm ~=0 )
 * where L=linear operator (Fourier matrix), D=regularization operator (identity). 
 * The above problem can be recast as:   
 *         min_m |G m -dd| 
 * where G = [L, eps D]^t and dd=[d,0]. 
 *-----------------------------------------------------------------------
 *
 * Copyright (c) 2021 Harbin Institute of Technology. All rights reserved.
 * Anothr: Pengliang Yang 
 * Email: ypl.2100@gmail.com
 * Homepage: https://yangpl.wordpress.com
 *-----------------------------------------------------------------------*/
#include "cstd.h"

void sf_matmult_init(float omega0, float dt_, float *omegas,int nfreq_,int nt_, int ntr_);
void sf_matmult_lop (bool adj, bool add, int nx, int ny, float *x, float *y);


void sf_cgstep( bool forget             /* restart flag */, 
		int nx                  /* model size */, 
		int ny                  /* data size */, 
		float * x        /* current model [nx] */,  
		const float * g  /* gradient [nx] */, 
		float * rr       /* data residual [ny] */,
		const float * gg /* conjugate gradient [ny] */) ;
/*< Step of Claerbout's conjugate-gradient iteration for complex operators. 
  The data residual is rr = A x - dat
  >*/

void sf_cgstep_close (void) ;
/*< Free allocated space. >*/ 

int nd, nm, nreg, niter, nrestart;
bool verb;
float tol, eps;
float *gr, *mm, *gm, *rr, *dd;

void adjsrc_freq2time_init(float omega0, float dt, float *omegas, int nfreq, int nt, int ntr,
			   float tol_, float eps_, int nrestart_, int niter_, int verb_)
{
  sf_matmult_init(omega0, dt, omegas, nfreq, nt, ntr);

  nd = 2*nfreq*ntr;
  nm = nt*ntr;
  nreg = nm;
  tol = tol_;/* tolerance for convergence */
  eps = eps_;/* regularization coefficient */
  nrestart = nrestart_; /* restart after every nrestart iterations */
  niter = niter_;
  verb = verb_?true:false;

  rr = alloc1float(nd+nreg); /* rr=[Lm-d, eps*Dm]^t*/
  gm = alloc1float(nm);/* gradient w.r.t. model mm */
  gr = alloc1float(nd+nreg);/* conjugate gradient w.r.t. residual rr */
}

void adjsrc_freq2time_close()
{
  free(rr);
  free(gm);
  free(gr);
}


void adjsrc_init(float *mm_, float *dd_)
/*< assign pointers to data (freq domain emf) and model (time domain emf)>*/
{
  /* here we avoid repeated allocation of memory */
  mm = mm_;
  dd = dd_;
}


void adjsrc_freq2time()
/*< perform many iterations of least-squares optimization >*/
{
  bool forget=false;

  int i, iter;
  float dpr0,dpg0,dpr,dpg,res_mod;

  for(i=0; i<nm; i++)   mm[i] = 0.;/* initialize model with 0 vector: m=0 */
  for(i=0; i<nd; i++)   rr[i] =- dd[i];/* rr_part1 = Lm-d=-d when m=0 */
  for(i=0; i<nreg; i++) rr[i+nd] = 0.; /* rr_part2 = eps*Dm=0 when m=0 */

  for(dpr0=0., i=0; i<nd; i++) dpr0 += rr[i]*conj(rr[i]);
  dpg0 = 1.;
  for(iter=0; iter<niter; iter++) {
    /* apply adjoint of operator G to residual rr: G^t rr=[L^t,eps*I]rr */
    sf_matmult_lop(true,  false, nm, nd, gm, rr);/* gm=Lt [rr_part1] */
    for(i=0; i<nreg; i++) gm[i] += eps*rr[i+nd]; // gm += Dt[rr_part2], D=Identity here

    /* apply forward operator G to gradient gm: gr=[L,eps*I]gm */
    sf_matmult_lop(false, false, nm, nd, gm, gr);// gr_part1=L[gm]
    for(i=0; i<nreg; i++) gr[i+nd] = eps*gm[i];/* gr_part2=eps*gm */

    if(forget && nrestart!=0) forget = (bool)(0==(iter+1)%nrestart);

    if (iter == 0) {
      for(dpg0=0.,i=0; i<nm; i++) dpg0 += gm[i]*conj(gm[i]);
      dpr = 1.;
      dpg = 1.;
    } else {
      for(dpr=0., i=0; i<nd; i++) dpr += rr[i]*conj(rr[i]);
      dpr /= dpr0;

      for(dpg=0., i=0; i<nm; i++) dpg += gm[i]*conj(gm[i]);
      dpg /= dpg0;
    }
    if (verb){
      for(res_mod=0., i=nd; i<nd+nreg; i++) res_mod += rr[i]*conj(rr[i]);
      printf("iteration %d res dat %e res mod %e grad %e\n", iter, dpr, res_mod, dpg);
    }
    if (dpr < tol || dpg < tol) {
      if (verb)	printf("convergence in %d iterations\n",iter+1);
      break;
    }

    /* Claerbout's CG: (mm, rr)=cgstep(mm, rr, gm, gr); */
    sf_cgstep(forget, nm, nd+nreg, mm, gm, rr, gr); 
    forget = false;
  }
  sf_cgstep_close();

  /* least-squares solution is now in vector mm */
}
